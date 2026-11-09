/* Build a clean, non-destructive GitHub-ready release folder for the insomnia TCC.
 * It never edits the working project. Re-run with --overwrite to rebuild output.
 */
const fs = require('fs');
const path = require('path');
const root = path.resolve(__dirname, '..');
const p = (...x) => path.join(root, ...x);
const dest = p('tcc_insomnia_gwas_github');
const overwrite = process.argv.includes('--overwrite');
if (fs.existsSync(dest) && !overwrite) throw new Error(`${dest} already exists. Re-run with --overwrite to rebuild it.`);
if (fs.existsSync(dest)) fs.rmSync(dest, { recursive: true, force: true });
fs.mkdirSync(dest, { recursive: true });
const cp = (source, target) => { if (!fs.existsSync(source)) { console.warn(`Missing, skipped: ${source}`); return; } fs.mkdirSync(path.dirname(target), {recursive:true}); fs.copyFileSync(source, target); };
const text = (target, value) => { fs.mkdirSync(path.dirname(target), {recursive:true}); fs.writeFileSync(target, value); };
const project = p('projeto_insonia_tcc');
const analytic = path.join(project, 'Resultados', 'analise_padronizada');

// Final executable scripts. Earlier exploratory scripts and the obsolete FinnGen R20 integration are intentionally excluded.
['deduplicate_finngen_r13_fuma_input.py','process_finngen_r13_to_fuma.py','process_opengwas_vcf_to_fuma.py','process_sumstats_to_fuma.py','validate_finngen_r13_fuma_input.py','validate_fuma_inputs.py','validate_opengwas_fuma_inputs.py'].forEach(f => cp(path.join(project,'pipeline_fuma',f), path.join(dest,'scripts_em_python',f)));
['17_compare_watanabe_supplement.py','18_compare_gwas_catalog_insomnia.py','20_prepare_watanabe_hcp_context_inputs.py','21_prepare_finngen_r13_corrected_enrichment_inputs.py','22_integrate_final_enrichment_results.py','23_integrate_finngen_r13_corrected_core16.py','24_prepare_gwas_catalog_global_inputs.py','25_integrate_string_results.js','26_assemble_github_release.js','27_build_tcc_figures_and_update_docx.py'].forEach(f => cp(p('scripts',f), path.join(dest,'scripts_em_python_e_node',f)));
['08_analyze_fuma_results.R','09_integrate_gene_evidence.R','10_compare_fuma_loci.R','11_compare_legacy_jansen_watanabe.R','12_integrate_opengwas_ukb_sensitivity.R','13_prepare_gene2func_inputs.R','14_prepare_final_gene_evidence_table.R','15_prepare_legacy_jansen_watanabe_sensitivity_inputs.R','16_integrate_magma_gene_sets.R','19_integrate_opengwas_ukb_gene_sets.R','21_plot_core16_evidence_heatmap.R'].forEach(f => cp(p('R',f), path.join(dest,'scripts_em_R',f)));

// Minimal reusable inputs and final results; no raw sumstats or full FUMA archives.
[
  '22_core16_final_evidence_table.tsv','15_cross_study_locus_clusters.tsv','16_pairwise_locus_overlap_summary.tsv','21_core16_opengwas_ukb_sensitivity_summary.tsv','30_core16_x_watanabe_supplement_comparison.tsv','32_gwas_catalog_insomnia_x_primary_fuma_loci.tsv','33_gwas_catalog_insomnia_locus_summary.tsv','47_registro_integrado_analises_funcionais.tsv','48_termos_funcionais_significativos_integrados.tsv','49_core16_x_finngen_r13_corrigido_magma.tsv','50_core16_finngen_r13_corrigido_validation_summary.tsv','51_string_network_summary.tsv','52_string_top_nodes.tsv','53_string_key_functional_terms.tsv','54_string_interpretation_notes.tsv'
].forEach(f => cp(path.join(analytic,f), path.join(dest,'results','tables',f)));
[
  ['enrichment_inputs','legacy_jansen_watanabe_28','legacy_jansen_watanabe_28_gene_symbols.txt'],
  ['enrichment_inputs','gwas_catalog_insomnia_global','gwas_catalog_insomnia_global_ensembl_ids.txt'],
  ['string_inputs','README_STRING.md']
].forEach(parts => cp(path.join(project,...parts), path.join(dest,'inputs',parts.slice(1).join('__'))));
const figures = [
  ['STRING_core16','STRING_network_core16.png'], ['STRING_core16','string_vector_graphic_core16.svg'],
  ['STRING_legacy28','STRING_network_legacy28.png'],
  ['STRING_global_gwas_catalog','STRING_network_global_gwas_catalog.png'], ['STRING_global_gwas_catalog','STRING_vector_graphic_network_global_gwas_catalog.svg'],
  ['STRING_global_gwas_catalog','enrichment_Process_sim0.8_graph_global_gwas_catalog.png'], ['STRING_global_gwas_catalog','enrichment_Process_sim0.8_graph_global_gwas_catalog.svg']
];
figures.forEach(([d,f]) => cp(path.join(project,'Resultados_STRING',d,f), path.join(dest,'results','figures',`${d}__${f}`)));
const rawString = [
  ['STRING_core16','string_interactions_core16.tsv'],['STRING_core16','string_node_degrees_core16.tsv'],['STRING_core16','enrichment.all_core16.tsv'],
  ['STRING_legacy28','STRING_interactions_legacy28.tsv'],['STRING_legacy28','STRING_node_degrees_legacy28.tsv'],['STRING_legacy28','enrichment.all_legacy28.tsv'],
  ['STRING_global_gwas_catalog','string_interactions.tsv'],['STRING_global_gwas_catalog','string_node_degrees_global_gwas_catalog.tsv'],['STRING_global_gwas_catalog','enrichment.all_global_gwas_catalog.tsv'],['STRING_global_gwas_catalog','enrichment.Process_global_gwas_catalog.tsv']
];
rawString.forEach(([d,f]) => cp(path.join(project,'Resultados_STRING',d,f), path.join(dest,'results','string_exports',`${d}__${f}`)));
['DOCUMENTACAO_COMPLETA_METODOLOGIA_E_SCRIPTS.md','DOCUMENTACAO_METODOLOGIA_ATE_GENE2FUNC.md','METODOLOGIA_POS_GENE2FUNC_E_PLANO_STRING.md','GUIA_GPROFILER_E_SHINYGO.md','GUIA_ESCRITA_TCC_METODOLOGIA_E_RESULTADOS.md'].forEach(f => cp(path.join(project,f), path.join(dest,'docs',f)));
cp(path.join(analytic,'DOCUMENTACAO_STRING_E_INTERPRETACAO.md'),path.join(dest,'docs','DOCUMENTACAO_STRING_E_INTERPRETACAO.md'));

text(path.join(dest,'.gitignore'), `# Raw or bulky source data: download from original repositories, never commit here\n*.gz\n*.zip\n*.vcf\n*.bcf\n*.RData\n*.Rhistory\nR_libs/\nsumstats_brutas/\nFUMA_job*/\nResultados_brutos/\n.DS_Store\nThumbs.db\n`);
text(path.join(dest,'data_availability','README.md'), `# Data availability\n\nThis repository contains code, parameters, small derived tables, and selected figures. It intentionally excludes raw GWAS summary statistics, full FUMA job archives, and large downloaded files. Obtain raw data from their primary repositories and follow each source's terms of use.\n\nThe final FinnGen validation uses the corrected R13 FUMA run (job 763972); earlier FinnGen outputs are not included as replication evidence.\n`);
text(path.join(dest,'README.md'), `# Arquitetura genética da insônia: pipeline reprodutível\n\nEste repositório contém a versão organizada do pipeline do TCC: preparação/QC de sumstats, FUMA SNP2GENE e MAGMA, integração multiestudo, enriquecimento funcional e redes STRING.\n\n## Resultado principal\n\nO conjunto **Core16** reúne genes protein-coding com significância MAGMA por Bonferroni em pelo menos duas unidades de evidência independentes e ao menos uma modalidade de mapeamento SNP2GENE. Ele está em \`results/tables/22_core16_final_evidence_table.tsv\`. Jansen e Watanabe foram tratados como uma única unidade na regra de consenso devido à sobreposição amostral.\n\n## Organização\n\n- \`scripts_python/\`: conversão, QC e validação de inputs FUMA.\n- \`scripts_python_and_node/\`: integrações complementares, incluindo GWAS Catalog, FinnGen corrigido e STRING.\n- \`scripts_R/\`: integração dos resultados FUMA/MAGMA.\n- \`inputs/\`: listas pequenas reutilizáveis para análises de sensibilidade e STRING.\n- \`results/tables/\`: tabelas finais para o manuscrito.\n- \`results/figures/\`: figuras STRING selecionadas.\n- \`results/string_exports/\`: exportações tabulares originais do STRING que sustentam a integração.\n- \`docs/\`: metodologia, interpretação e guia de escrita.\n- \`data_availability/\`: política para dados brutos não versionados.\n\n## Como interpretar as redes\n\nAs três redes correspondem a três níveis de evidência: Core16 (principal), Legacy28 (sensibilidade histórica) e GWAS Catalog global (camada descritiva). A última revela convergência para sinalização sináptica, transmissão trans-sináptica, desenvolvimento do sistema nervoso e adesão celular, mas não redefine o Core16 nem demonstra causalidade. Veja \`results/tables/51_string_network_summary.tsv\` e \`docs/DOCUMENTACAO_STRING_E_INTERPRETACAO.md\`.\n\n## Dados não incluídos\n\nNão suba estatísticas-resumo brutas, dados individuais, arquivos VCF/BCF, bibliotecas R, ZIPs completos do FUMA ou screenshots de interface. Eles permanecem no diretório de trabalho local e devem ser recuperados nas fontes oficiais, conforme licenças e versões descritas em \`docs/\`.\n\n## Antes de tornar público\n\nRevise autoria, afiliações, referências, licença escolhida pelo grupo e se as políticas de cada repositório permitem redistribuir cada arquivo derivado.\n`);
const readmePath = path.join(dest, 'README.md');
fs.writeFileSync(readmePath, fs.readFileSync(readmePath, 'utf8')
  .replaceAll('scripts_python/', 'scripts_em_python/')
  .replaceAll('scripts_python_and_node/', 'scripts_em_python_e_node/')
  .replaceAll('scripts_R/', 'scripts_em_R/'));
console.log(`Created clean GitHub release at ${dest}`);
