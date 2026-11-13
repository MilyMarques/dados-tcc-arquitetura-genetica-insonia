# Guia de escrita: Metodologia e Resultados

Este guia organiza os resultados finais do TCC sem alterar a hierarquia de evidências: **Core16** é o resultado principal; **Legacy28**, **GWAS Catalog global**, HCP289, UKB OpenGWAS e FinnGen R13 são camadas de sensibilidade, contexto ou validação externa.

## Metodologia

### 1. Desenho do estudo e fontes

Descreva um estudo de bioinformática baseado em estatísticas-resumo públicas de GWAS de insônia. Liste os quatro estudos primários (Hammerschlag, Lane, Jansen e Watanabe), suas coortes, ancestralidade, definição fenotípica, tamanho amostral e fonte dos dados. A tabela-base é `Resultados/analise_padronizada/01_study_summary.tsv`.

Explique separadamente as camadas não primárias:

- UKB OpenGWAS: análises de sensibilidade (`21_core16_opengwas_ukb_sensitivity_summary.tsv`);
- FinnGen R13: validação externa do endpoint clínico finlandês (`50_core16_finngen_r13_corrigido_validation_summary.tsv`);
- GWAS Catalog: camada descritiva de genes mapeados em associações publicadas, não uma meta-análise nem evidência causal;
- Legacy28: análise histórica de sensibilidade, derivada das primeiras rodadas Jansen–Watanabe;
- HCP289 de Watanabe: contexto funcional publicado, sem independência em relação ao estudo de origem.

### 2. Pré-processamento, QC e submissão ao FUMA

Informe que as estatísticas foram padronizadas, auditadas e validadas antes da geração do input compactado para FUMA. Declare a construção genômica usada em cada input; o FinnGen R13 foi corrigido para a construção exigida pelo FUMA antes da reanálise. Nomeie os scripts de `scripts_em_python/` como a rota final de processamento. Os scripts R `00_config.R` a `07_run_study.R`, presentes apenas no diretório histórico de trabalho, não são a rota final e não integram este release.

### 3. FUMA SNP2GENE e MAGMA

Reporte os parâmetros efetivamente fixados entre os quatro estudos principais: painel 1000G Phase 3 EUR; SNPs independentes significativos com r² <0,6; lead SNPs com r² <0,1; limiar de associação P<5×10⁻⁸; SNPs candidatos P<0,05; janela de locus 250 kb; MHC excluído do SNP2GENE; e mapeamento posicional, eQTL e interação cromatínica conforme os jobs. Para os outputs finais integrados, registre FUMA v1.8.2, MAGMA v1.08 e Ensembl v110_userGen; as configurações históricas v1.8.3/v1.10/v110 de Jansen e Watanabe devem ser mencionadas, se necessário, apenas como execuções anteriores preservadas no arquivo. Para MAGMA, diferencie claramente:

- análise **gene-based**, que produz `magma.genes.out` e fornece a significância Bonferroni por gene;
- análise **gene-set**, que testa conjuntos/vias e não selecionou os genes do Core16.

Os resultados gene-based e evidências FUMA foram integrados nas tabelas `04_mapped_genes_and_evidence.tsv`, `05_magma_all_genes.tsv`, `06_magma_bonferroni_genes.tsv`, `09_gene_x_evidence_x_study.tsv` e `10_gene_evidence_summary.tsv`.

### 4. Definição do conjunto principal

Declare previamente a regra operacional: gene *protein-coding*, significativo por Bonferroni no MAGMA em pelo menos duas unidades analíticas e com pelo menos uma modalidade SNP2GENE (posicional, eQTL ou interação cromatínica). Jansen e Watanabe foram tratados como uma única unidade UK Biobank devido à sobreposição de amostras. A lista final está em `22_core16_final_evidence_table.tsv`.

### 5. Comparação de loci, enriquecimento e redes

Explique que loci foram agrupados por proximidade/partilha entre estudos (`15_cross_study_locus_clusters.tsv`, `16_pairwise_locus_overlap_summary.tsv`). Para GENE2FUNC, g:Profiler, ShinyGO e STRING, deixe explícito o papel de cada lista e o controle de múltiplos testes usado por cada plataforma. Diferencie:

- Core16: interpretação de uma seleção replicada e conservadora;
- Legacy28: sensibilidade histórica;
- global do GWAS Catalog: caracterização funcional descritiva de associações publicadas.

No STRING: *Homo sapiens*, rede completa, sem nós adicionais, escore >=0,700 e FDR 1%. Descreva a rede como análise de interação conhecida/predita, não como prova experimental de interação ou causalidade.

## Resultados: sequência sugerida

1. **Estudos e QC.** Tabela de estudos e fluxograma de inclusão/submissão.
2. **Mapeamento e genes MAGMA.** Resumo de genes Bonferroni por estudo e quadro de evidências SNP2GENE.
3. **Loci compartilhados.** Tabela resumida de clusters/loci entre estudos.
4. **Core16.** Tabela principal com símbolo, Ensembl, estudos MAGMA e modalidades de mapeamento. Apresente-a como resultado central.
5. **Sensibilidade/validação.** UKB OpenGWAS e FinnGen R13 em uma tabela concisa; FinnGen corrigido não replicou estritamente o Core16 e deve ser reportado como resultado negativo informativo.
6. **Função e vias.** Explique que listas pequenas (Core16/Legacy28) não tiveram enriquecimento robusto replicado; a camada global do GWAS Catalog mostrou convergência descritiva para sinapse, transmissão trans-sináptica, desenvolvimento neural e adesão celular. Use `47_registro_integrado_analises_funcionais.tsv` e `48_termos_funcionais_significativos_integrados.tsv`.
7. **STRING.** Apresente as três redes como um funil: global -> Legacy28 -> Core16. A tabela comparativa está em `51_string_network_summary.tsv`; os termos do STRING em `53_string_key_functional_terms.tsv`.

## Figuras e tabelas recomendadas para o corpo do TCC

### Tabelas principais

1. Características dos GWAS primários: `01_study_summary.tsv`.
2. Core16 e suas evidências: `22_core16_final_evidence_table.tsv`.
3. Resumo de loci compartilhados: `16_pairwise_locus_overlap_summary.tsv` (versão condensada).
4. Sensibilidade/validação: combinar `21_core16_opengwas_ukb_sensitivity_summary.tsv` e `50_core16_finngen_r13_corrigido_validation_summary.tsv`.
5. Integração funcional/STRING: `51_string_network_summary.tsv` e uma seleção não redundante de `53_string_key_functional_terms.tsv`.

### Figuras principais

1. Fluxograma do pipeline: fontes -> QC -> FUMA/MAGMA -> Core16 -> camadas de sensibilidade/contexto -> enriquecimento/STRING.
2. Esquema/heatmap de evidências do Core16, derivado de `22_core16_final_evidence_table.tsv`.
3. Diagrama de loci compartilhados ou UpSet plot, derivado de `15_cross_study_locus_clusters.tsv`.
4. Painel STRING com as três camadas: `Resultados_STRING/STRING_global_gwas_catalog/STRING_vector_graphic_network_global_gwas_catalog.svg`, `Resultados_STRING/STRING_legacy28/STRING_network_legacy28.png`, e `Resultados_STRING/STRING_core16/string_vector_graphic_core16.svg`.
5. Grafo de processos STRING global: `Resultados_STRING/STRING_global_gwas_catalog/enrichment_Process_sim0.8_graph_global_gwas_catalog.svg`.

Deixe screenshots de interface e tabelas extensas em apêndices ou material suplementar. Não use “Insomnia measurement”/“sleep” do STRING como descoberta funcional: são anotações circulares porque a lista foi montada com associações de insônia.

## Apêndices e reprodutibilidade

Inclua ou deposite no GitHub: todos os scripts finais, parâmetros FUMA, manifestos/README, listas de genes, tabelas TSV finais e figuras vetoriais. Não suba estatísticas-resumo brutas, grandes ZIPs de jobs FUMA ou dados licenciados; cite a fonte e forneça instruções de download.
