/*
 * Integrates the three manually exported STRING networks for the insomnia TCC.
 * Input: project_insonia_tcc/Resultados_STRING/{STRING_core16,STRING_legacy28,STRING_global_gwas_catalog}
 * Output: project_insonia_tcc/Resultados/analise_padronizada/51-54_*.tsv and DOCUMENTACAO_STRING_E_INTERPRETACAO.md
 * Requires Node.js >=18; no external packages.
 */
const fs = require('fs');
const path = require('path');

const root = path.resolve(__dirname, '..');
const project = path.join(root, 'projeto_insonia_tcc');
const stringRoot = path.join(project, 'Resultados_STRING');
const out = path.join(project, 'Resultados', 'analise_padronizada');
fs.mkdirSync(out, { recursive: true });

const networks = [
  { id: 'core16', label: 'Core16 (principal)', dir: 'STRING_core16', inputs: 16, mapped: 16, ppi: '0.0573', role: 'lista principal, priorizada por evidência multiestudo' },
  { id: 'legacy28', label: 'Legacy28 (sensibilidade histórica)', dir: 'STRING_legacy28', inputs: 28, mapped: 27, ppi: '0.0153', role: 'sensibilidade histórica; um gene não codificante não foi mapeado pelo STRING' },
  { id: 'global_gwas_catalog', label: 'GWAS Catalog global (descritiva)', dir: 'STRING_global_gwas_catalog', inputs: 786, mapped: 643, ppi: 'NA (estatística não exportada)', role: 'camada descritiva de genes mapeados em associações publicadas de insônia' }
];

const fileFor = (n, stem) => {
  const d = path.join(stringRoot, n.dir);
  const exact = fs.readdirSync(d).find(x => x.toLowerCase() === stem.toLowerCase());
  if (exact) return path.join(d, exact);
  return null;
};
const readTsv = f => {
  const lines = fs.readFileSync(f, 'utf8').replace(/^\uFEFF/, '').trim().split(/\r?\n/);
  return lines.map(x => x.split('\t'));
};
const q = x => String(x ?? '').replace(/[\t\r\n]/g, ' ').trim();
const writeTsv = (name, header, rows) => fs.writeFileSync(path.join(out, name), [header, ...rows].map(r => r.map(q).join('\t')).join('\n') + '\n');

const summary = [];
const nodeRows = [];
const termRows = [];
for (const n of networks) {
  const d = path.join(stringRoot, n.dir);
  const degreeFile = fs.readdirSync(d).find(x => /^string_node_degrees.*\.tsv$/i.test(x));
  const interactionFile = fs.readdirSync(d).find(x => /^string_interactions.*\.tsv$/i.test(x) && !/psimi/i.test(x));
  const degrees = readTsv(path.join(d, degreeFile));
  const degreeHeader = degrees[0].map(x => x.replace(/^#/, '').toLowerCase());
  const nodeIdx = degreeHeader.findIndex(x => x === 'node');
  const degIdx = degreeHeader.findIndex(x => x === 'node_degree');
  const degreeValues = degrees.slice(1).map(r => ({ node: r[nodeIdx], degree: Number(r[degIdx]) || 0 })).sort((a,b) => b.degree - a.degree || a.node.localeCompare(b.node));
  const interactionRows = readTsv(path.join(d, interactionFile));
  const h = interactionRows[0].map(x => x.replace(/^#/, '').toLowerCase());
  const a = h.indexOf('node1'), b = h.indexOf('node2'), score = h.indexOf('combined_score');
  const edges = new Map();
  for (const r of interactionRows.slice(1)) {
    const x = r[a], y = r[b], s = Number(r[score]);
    if (!x || !y || s < 0.7) continue;
    const k = [x,y].sort().join('::');
    if (!edges.has(k)) edges.set(k, s);
  }
  const sortedEdges = [...edges.entries()].sort((x,y) => y[1]-x[1]);
  const edgeText = sortedEdges.slice(0, 20).map(([k,v]) => `${k.replace('::','–')} (${v.toFixed(3)})`).join('; ') + (sortedEdges.length > 20 ? `; … (+${sortedEdges.length - 20} arestas no arquivo bruto)` : '');
  const edgeCount = edges.size;
  const nonIsolated = degreeValues.filter(x => x.degree > 0).length;
  const meanDegree = n.mapped ? (2*edgeCount/n.mapped) : 0;
  summary.push([n.id, n.label, n.role, n.inputs, n.mapped, n.inputs-n.mapped, edgeCount, nonIsolated, n.mapped-nonIsolated, meanDegree.toFixed(3), n.ppi, edgeText || 'nenhuma aresta >= 0.700']);
  degreeValues.slice(0, 20).forEach((x, rank) => nodeRows.push([n.id, n.label, rank+1, x.node, x.degree]));

  const enrichName = n.id === 'global_gwas_catalog' ? 'enrichment.all_global_gwas_catalog.tsv' : `enrichment.all_${n.id}.tsv`;
  const enrichFile = fileFor(n, enrichName);
  if (enrichFile) {
    const tab = readTsv(enrichFile);
    const eh = tab[0].map(x => x.replace(/^#/, '').toLowerCase());
    const find = x => eh.indexOf(x);
    const cat = find('category'), termId = find('term id'), desc = find('term description'), obs = find('observed gene count'), fdr = find('false discovery rate'), labels = find('matching proteins in your network (labels)');
    const reportableCategories = new Set(['GO Process','GO Component','Reactome','KEGG','Local Network Cluster']);
    const accepted = tab.slice(1).map(r => ({ r, f: Number(r[fdr]) })).filter(x => reportableCategories.has(x.r[cat]) && Number.isFinite(x.f) && x.f <= 0.01);
    for (const category of reportableCategories) {
      accepted.filter(x => x.r[cat] === category).sort((x,y)=>x.f-y.f).slice(0,5).forEach(x => {
        const r=x.r; termRows.push([n.id, r[cat] || 'NA', r[termId] || 'NA', r[desc] || 'NA', r[obs] || 'NA', r[fdr] || 'NA', r[labels] || '']);
      });
    }
  }
}

writeTsv('51_string_network_summary.tsv', ['network_id','network_label','analytical_role','input_genes','string_mapped_proteins','unmapped_input_genes','undirected_edges_score_ge_0.700','connected_nodes','isolated_nodes','mean_degree','STRING_PPI_enrichment_p','edges_score_ge_0.700'], summary);
writeTsv('52_string_top_nodes.tsv', ['network_id','network_label','rank_by_degree','gene_symbol','node_degree'], nodeRows);
writeTsv('53_string_key_functional_terms.tsv', ['network_id','category','term_id','term_description','observed_gene_count','STRING_FDR','matching_gene_symbols'], termRows);
writeTsv('54_string_interpretation_notes.tsv', ['network_id','interpretive_note'], [
  ['core16','A rede principal é muito esparsa: somente HEXIM1–HEXIM2 permaneceu com score >=0.700. PPI enrichment p=0.0573, portanto não há evidência formal de conectividade acima do acaso a 0.05.'],
  ['legacy28','A sensibilidade histórica mostrou duas arestas (HEXIM1–HEXIM2 e BTBD9–MEIS1) e PPI enrichment p=0.0153. Isso indica conectividade acima da expectativa, mas não um módulo denso; é uma análise de sensibilidade, não substitui o Core16.'],
  ['global_gwas_catalog','A rede global é descritiva e apresenta convergência para sinalização sináptica, transmissão trans-sináptica, desenvolvimento do sistema nervoso e adesão celular. Ela não demonstra causalidade, porque os genes são mapeamentos de associações heterogêneas do GWAS Catalog.'],
  ['all','As anotações de fenótipo “insônia/sleep” no STRING são circulares para listas selecionadas por associações de insônia e não devem ser tratadas como descoberta funcional.']
]);

const md = `# STRING: integração e interpretação\n\n## Configuração comum\n\nAs três redes foram geradas manualmente no STRING para *Homo sapiens*, tipo **full STRING network**, sem nós adicionais, com escore de confiança alto (>=0,700) e controle de FDR estrito (1%). O arquivo \`51_string_network_summary.tsv\` é a síntese quantitativa; \`52_string_top_nodes.tsv\` lista os nós de maior grau; \`53_string_key_functional_terms.tsv\` preserva os termos com FDR <=0,01.\n\n## Leitura por camada\n\n- **Core16:** 16 proteínas mapeadas, 1 aresta (HEXIM1–HEXIM2; score 0,995), PPI enrichment *p*=0,0573. A figura mostra heterogeneidade funcional entre os genes priorizados; não sustenta um único módulo proteico. É a figura de rede principal porque os genes foram definidos com a regra mais rigorosa.\n- **Legacy28:** 27/28 genes mapeados (o transcrito não codificante RP11-463D19.2 não possui proteína humana correspondente), 2 arestas e PPI enrichment *p*=0,0153. A conexão adicional é BTBD9–MEIS1 (score 0,765). Trate como análise de sensibilidade histórica, não como validação independente.\n- **Global GWAS Catalog:** 643/786 entradas foram mapeadas a proteínas, com 359 arestas. Esta escala revela a arquitetura funcional descritiva de genes associados à insônia em estudos publicados: sinalização sináptica, transmissão trans-sináptica, desenvolvimento do sistema nervoso e adesão celular. A imagem é visualmente informativa, mas não pode redefinir o Core16.\n\n## Como reportar\n\nUse o painel das três redes para explicar o funil analítico: ampla evidência publicada (global) -> seleção histórica (Legacy28) -> conjunto replicado e priorizado (Core16). Não conclua que a rede global prova mecanismo causal ou que o PPI da Legacy28 represente uma via única.\n`;
fs.writeFileSync(path.join(out, 'DOCUMENTACAO_STRING_E_INTERPRETACAO.md'), md);
console.log(`Wrote STRING integration tables to ${out}`);
