# TCC — Caracterização da Arquitetura Genética da Insônia

> **Nota para este release GitHub.** Este documento preserva nomes de caminhos do diretório de trabalho original (`R/`, `scripts/`, `Resultados/`). No release, os equivalentes estão em `scripts_em_R/`, `scripts_em_python_e_node/` e `results/tables/`. Estatísticas-resumo brutas e ZIPs FUMA não são redistribuídos; por isso, o release documenta e audita o fluxo, mas não permite execução integral sem os downloads originais.

## Escopo e regra de interpretação

Este projeto integra estatísticas-resumo de GWAS de insônia, executa mapeamento SNP–gene e análise MAGMA no FUMA, e realiza interpretação funcional posterior. A cadeia principal de evidência inclui Hammerschlag 2017, Lane 2019 (fenótipo extremo), Jansen 2019 (componente UK Biobank público) e Watanabe 2022 (componente UK Biobank público).

Jansen e Watanabe não são unidades inteiramente independentes: ambos representam dados públicos UK Biobank, enquanto as meta-análises completas publicadas incorporaram dados 23andMe que não estavam disponíveis neste projeto. Assim, para o consenso, Jansen e Watanabe contam como uma unidade de evidência `ukb_public`. Os dois OpenGWAS UKB são sensibilidade técnica, e FinnGen R13 é validação externa em um endpoint clínico finlandês.

O núcleo principal (`core16`) não foi escolhido por enriquecimento funcional. Foi definido antes dessas análises como 16 genes *protein-coding* que satisfazem simultaneamente: (i) significância MAGMA com Bonferroni em pelo menos duas unidades analíticas; e (ii) ao menos uma modalidade positiva de mapeamento SNP2GENE (posicional, eQTL cerebral ou interação cromatínica cortical). Jansen e Watanabe compõem uma única unidade UK Biobank, portanto essa regra não pressupõe independência completa de participantes. A tabela-fonte é `Resultados/analise_padronizada/22_core16_final_evidence_table.tsv`.

## Hierarquia relevante de pastas

```text
projeto_codex/
├── R/                                  # scripts R de integração/análise
├── scripts/                            # scripts Python de extensão e integração final
├── projeto_insonia_tcc/
│   ├── sumstats_brutas/                # arquivos baixados, preservados sem alteração
│   ├── fuma_inputs_final/              # inputs FUMA principais, compactados
│   ├── fuma_inputs_opengwas/           # inputs de VCF OpenGWAS
│   ├── fuma_inputs_external/           # input FinnGen R13 e auditorias
│   ├── Resultados/                     # ZIPs FUMA dos quatro estudos e análises padronizadas
│   │   └── analise_padronizada/        # tabelas reprodutíveis de integração
│   ├── Resultados_finn/                # jobs FinnGen; o ZIP corrigido é o job 763972
│   ├── Resultados_opengwas/            # jobs UKB OpenGWAS
│   ├── Resultados_GENE2FUNC/           # ZIPs de GENE2FUNC
│   ├── Resultados_ShinyGO/             # tabelas/figuras ShinyGO
│   ├── Resultados_gProfiler/           # tabelas/figuras g:Profiler
│   ├── enrichment_inputs/              # query e background de cada análise funcional
│   └── string_inputs/                  # lista para a futura rede STRING
└── scripts_validacao/                  # verificações históricas de entradas/resultados
```

## Metodologia executada

### 1. Inventário e obtenção dos estudos

Foram registrados fenótipo, ancestralidade, tamanho amostral, fonte, versão do arquivo e build genômico. Estatísticas-resumo foram obtidas de repositórios públicos dos autores/CNCR/GWAS Catalog, OpenGWAS e FinnGen. O arquivo FinnGen R13 utilizado foi `sumstats_brutas/finngen_r13/finngen_R13_F5_INSOMNIA.gz`, endpoint clínico F5_INSOMNIA.

O arquivo histórico OpenGWAS `finn-b-F5-INSOMNIA.vcf` falhou no SNP2GENE por ausência de SNPs candidatos significativos e não é a fonte usada para a análise externa final. A análise final FinnGen usa o release R13 e o job FUMA corrigido.

### 2. Processamento e controle de qualidade das sumstats

O pipeline Python padroniza cada fonte a um esquema mínimo para FUMA: `SNP`, `CHR`, `BP`, `A1`, `A2`, `P` e `N`; quando disponíveis, `BETA`, `SE`, `EAF`, `INFO` e metadados são retidos nos arquivos de auditoria.

As verificações aplicadas incluem:

- cabeçalho, delimitador e compactação gzip válidos;
- cromossomos permitidos (autossomos 1–22) e posição inteira válida;
- P em (0,1], ausência de valores obrigatórios e duplicatas;
- alelos em formato A/C/G/T, sem valores ausentes;
- identificação do build e, no FinnGen R13, conversão GRCh38 → GRCh37/hg19 antes do FUMA;
- relatório JSON de contagens antes/depois e validação formal do `.FUMA.txt.gz`.

Os arquivos FUMA não são “estatísticas completas”: são versões padronizadas para submissão ao SNP2GENE. As estatísticas brutas e os JSONs de QC permanecem preservados para rastreabilidade.

### 3. SNP2GENE e MAGMA no FUMA

Cada input foi submetido ao FUMA SNP2GENE em hg19/GRCh37, com painel 1000G Phase 3 EUR para estudos de ancestralidade europeia. Foram adotados P de lead SNP = `5e-8`, P de SNP candidato = `0,05`, r² = `0,6` para SNPs independentes significativos, r² = `0,1` para lead SNPs e distância máxima entre blocos = 250 kb.

Foram habilitados mapeamento posicional, eQTL cerebral GTEx v8 e interação cromatínica de córtex adulto/fetal de Giusti-Rodríguez et al. O MHC foi excluído de anotações SNP2GENE, mantendo o procedimento consistente entre estudos. MAGMA foi executado pelo FUMA com janela de 10 kb. GTEx v8 e BrainSpan foram usados posteriormente para contextualização no GENE2FUNC, não como análise formal MAGMA gene-property.

Para genes MAGMA, a significância primária por estudo foi `P < 0,05 / número de genes testados`; FDR Benjamini–Hochberg é apresentado como suporte secundário. Para conjuntos MAGMA e GENE2FUNC, foi usado FDR BH < 0,05; para g:Profiler, g:SCS < 0,05; e para ShinyGO, FDR < 0,05.

### 4. Integração dos quatro estudos principais

Os ZIPs do FUMA foram lidos sem modificar os originais. Foram extraídos `leadSNPs.txt`, `GenomicRiskLoci.txt`, `genes.txt`, `eqtl.txt`, `ci.txt`, `magma.genes.out` e `magma.gsa.out`. A análise gerou:

- comparações de loci por interseção de intervalos FUMA no mesmo cromossomo;
- matriz gene × evidência × estudo;
- comparação de genes MAGMA e definição auditável do core16;
- comparação histórica entre a lista prévia Jansen–Watanabe e a análise atual;
- análise MAGMA gene-set genome-wide, mantida conceitualmente separada de enrichment por lista de genes.

As tabelas `09_gene_x_evidence_x_study.tsv`, `10_gene_evidence_summary.tsv`, `15_cross_study_locus_clusters.tsv`, `16_pairwise_locus_overlap_summary.tsv` e `22_core16_final_evidence_table.tsv` são as fontes principais para Resultados.

### 5. Sensibilidades e validações externas

Os jobs `ukb-a-13` e `ukb-b-3957` foram usados como sensibilidade UKB; por compartilharem a coorte, não modificam o core. A tabela `21_core16_opengwas_ukb_sensitivity_summary.tsv` descreve sua concordância.

FinnGen R13 foi submetido ao FUMA após correção do input. O resultado final é `Resultados_finn/FUMA_job763972_finnr13_corrigido-full.zip`; o ZIP `FUMA_job763398_finn_r13.zip` é anterior à correção e deve ser tratado como arquivo histórico, não como fonte final. A query funcional contém 21 genes abaixo do limiar MAGMA de Bonferroni de 0,05/18.802 = 2,659 × 10⁻⁶; o background funcional contém 18.225 genes com P MAGMA finito e símbolo utilizável. Na comparação corrigida, nenhum dos 16 genes do núcleo alcançou Bonferroni no FinnGen; isso é um resultado de ausência de replicação estrita, plausivelmente compatível com diferenças de fenótipo (insônia clínica), população e arquitetura genética. Não deve ser convertido em evidência de confirmação externa.

Watanabe HCP289 é uma análise de contexto da literatura: os 289 genes de alta confiança foram extraídos da tabela suplementar publicada e comparados ao background de genes mapeados nos loci de Watanabe. Como essa publicação contém UKB e 23andMe, ela amplia a interpretação funcional, mas não constitui validação independente nem redefine o core16.

### 6. Enriquecimento funcional

Foram analisadas quatro listas: core16 (principal), legacy28 Jansen–Watanabe (sensibilidade histórica), HCP289 de Watanabe (contexto publicado) e FinnGen R13 corrigido (validação externa). Para cada uma, a query foi comparada a um background compatível com o mecanismo de seleção, nunca ao conjunto indiscriminado de todos os genes conhecidos.

Os resultados finalizados estão registrados em `Resultados/analise_padronizada/47_registro_integrado_analises_funcionais.tsv`; os termos FDR-significativos estão em `48_termos_funcionais_significativos_integrados.tsv`.

Resultado central: o HCP289 apresentou sinais funcionais FDR-significativos no GENE2FUNC e ShinyGO, incluindo desenvolvimento/organização neuronal, assinaturas de neurônios do mesencéfalo/GABAérgicos, interações L1CAM, junções comunicantes e sinalização celular. Termos nomeados por doenças neurodegenerativas são interpretados apenas como sobreposição de genes/processos compartilhados, não como evidência de que insônia seja uma doença neurodegenerativa.

O g:Profiler não encontrou termos g:SCS-significativos nas listas. As exportações feitas com `All results` são nominais/exploratórias e não entram como resultado confirmado.

## Catálogo de scripts R

| Script | Finalidade e status |
|---|---|
| `R/00_config.R` | Configuração do pipeline R, builds, limiares, diretórios e metadados dos estudos. Legado; o processamento R completo foi descontinuado por desempenho. |
| `R/01_utils.R` | Funções de log, criação de diretórios, validação e estruturas vazias. Legado. |
| `R/02_read_audit.R` | Leitura/auditoria de tabelas e VCF em R. Legado. |
| `R/03_standardize.R` | Conversão de colunas para o padrão interno R. Legado. |
| `R/04_qc.R` | Filtros de QC em R (P, coordenadas, alelos, INFO/MAF e duplicatas). Legado. |
| `R/05_validate.R` | Métricas de validação, lambda GC e coerência beta/SE/P. Legado. |
| `R/06_export_fuma.R` | Exportação gzip para FUMA a partir do objeto R. Legado. |
| `R/07_run_study.R` | Orquestrador dos módulos R 02–06. Legado; não é a rota final de processamento. |
| `R/08_analyze_fuma_results.R` | Extrai e padroniza resultados dos quatro ZIPs FUMA principais; gera arquivos 01–07. Fonte final. |
| `R/09_integrate_gene_evidence.R` | Constrói matriz gene × estudo × modalidades de evidência; gera arquivos 09–12. Fonte final. |
| `R/10_compare_fuma_loci.R` | Compara intervalos FUMA por sobreposição; gera arquivos 13–16. Fonte final. |
| `R/11_compare_legacy_jansen_watanabe.R` | Documenta diferença entre resultados preliminares e análise atual; gera 17–18. Fonte final para sensibilidade histórica. |
| `R/12_integrate_opengwas_ukb_sensitivity.R` | Integra os dois jobs UKB OpenGWAS como sensibilidade; gera 19–21. Fonte final, mas não altera core. |
| `R/13_prepare_gene2func_inputs.R` | Prepara core16 e background MAGMA dos quatro estudos principais. Fonte final para GENE2FUNC core. |
| `R/14_prepare_final_gene_evidence_table.R` | Gera a tabela final core16, incorpora GENE2FUNC e prepara arquivos 22–24. Fonte final. |
| `R/15_prepare_legacy_jansen_watanabe_sensitivity_inputs.R` | Produz lista/background legacy28. Fonte final para sensibilidade. |
| `R/16_integrate_magma_gene_sets.R` | Integra `magma.gsa.out` dos quatro estudos; gera 25–29. Fonte final para gene-set genome-wide. |
| `R/19_integrate_opengwas_ukb_gene_sets.R` | Integra gene-set MAGMA dos dois UKB OpenGWAS; gera 35–39. Sensibilidade. |
| `R/20_integrate_finngen_r13_external_validation.R` | Script histórico apontando para `FUMA_job763398_finn_r13.zip` (pré-correção). Não usar para resultados finais; foi substituído pelos scripts Python 21 e 23. |
| `preparar_fuma.R`, `processando_watanabe.R`, `run_pipeline.R`, `analise_pos_fuma.R`, `testando.R` | Tentativas iniciais/diagnósticas. Arquivos preservados, mas não são parte do pipeline final. |

## Catálogo de scripts Python

| Script | Finalidade e status |
|---|---|
| `projeto_insonia_tcc/pipeline_fuma/process_sumstats_to_fuma.py` | Pipeline principal configurável para estatísticas-resumo tabulares: leitura em blocos, padronização, QC e geração `.FUMA.txt.gz`. Fonte final para os quatro estudos principais. |
| `projeto_insonia_tcc/pipeline_fuma/validate_fuma_inputs.py` | Valida headers, campos obrigatórios, P, cromossomos, coordenadas, alelos, duplicatas e compactação dos inputs principais. Fonte final. |
| `projeto_insonia_tcc/pipeline_fuma/process_opengwas_vcf_to_fuma.py` | Converte VCF OpenGWAS para input FUMA após interpretar campos FORMAT/INFO. Usado nos UKB OpenGWAS. |
| `projeto_insonia_tcc/pipeline_fuma/validate_opengwas_fuma_inputs.py` | Validador específico das saídas VCF OpenGWAS. |
| `projeto_insonia_tcc/pipeline_fuma/process_finngen_r13_to_fuma.py` | Processa as sumstats R13, incluindo preparação para liftOver. Fonte de processamento FinnGen. |
| `projeto_insonia_tcc/pipeline_fuma/deduplicate_finngen_r13_fuma_input.py` | Remove duplicatas após conversão de build e cria input FinnGen corrigido. Fonte final do input que gerou o job 763972. |
| `projeto_insonia_tcc/pipeline_fuma/validate_finngen_r13_fuma_input.py` | Valida especificamente o input FinnGen antes da submissão. |
| `projeto_insonia_tcc/pipeline_jansen/*` e `pipeline_watanabe/*` | Scripts iniciais específicos para Jansen/Watanabe e visualização; preservados para auditoria, substituídos pelo pipeline configurável. |
| `projeto_insonia_tcc/pipeline_summary_stats_fuma_magma_padronizado.py` | Versão consolidada inicial do processamento de sumstats; referência histórica. |
| `projeto_insonia_tcc/scripts_validacao/01–04_*.py` | Comparações históricas de inputs, top SNPs, loci e genes MAGMA Jansen/Watanabe. Auditoria complementar. |
| `scripts/17_compare_watanabe_supplement.py` | Extrai/contrasta dados suplementares de Watanabe com o core16. |
| `scripts/18_compare_gwas_catalog_insomnia.py` | Produz comparação externa entre associações do GWAS Catalog e loci FUMA. |
| `scripts/20_prepare_watanabe_hcp_context_inputs.py` | Gera query HCP289 e background contextual de Watanabe. Fonte final do HCP289. |
| `scripts/21_prepare_finngen_r13_corrected_enrichment_inputs.py` | Gera query de 21 genes e background de 18.225 genes diretamente do ZIP FinnGen corrigido. Fonte final. |
| `scripts/22_integrate_final_enrichment_results.py` | Gera as tabelas 47–48 de integração final de enrichment. Fonte final. |
| `scripts/23_integrate_finngen_r13_corrected_core16.py` | Compara core16 ao `magma.genes.out` do ZIP FinnGen corrigido e produz tabelas 49–50. Fonte final da validação FinnGen. |
| `projeto_insonia_tcc/testes_view.py` | Utilitário visual/diagnóstico; sem papel inferencial. |

## Arquivos a apresentar ao orientador

- `Resultados/analise_padronizada/22_core16_final_evidence_table.tsv` — genes núcleo, evidências e sensibilidade UKB.
- `Resultados/analise_padronizada/15_cross_study_locus_clusters.tsv` e `16_pairwise_locus_overlap_summary.tsv` — comparação formal dos loci.
- `Resultados/analise_padronizada/49_core16_x_finngen_r13_corrigido_magma.tsv` e `50_core16_finngen_r13_corrigido_validation_summary.tsv` — validação FinnGen correta; substituem os arquivos 43–46, que são históricos pré-correção.
- `Resultados/analise_padronizada/47_registro_integrado_analises_funcionais.tsv` — quadro geral de todas as rodadas funcionais.
- `Resultados/analise_padronizada/48_termos_funcionais_significativos_integrados.tsv` — termos significativos, genes contribuidores e limites interpretativos.
- `Resultados_GENE2FUNC/`, `Resultados_ShinyGO/` e `Resultados_gProfiler/` — resultados brutos e figuras exportadas por plataforma.

## Planejamento anterior (etapa concluída)

A rede STRING principal deveria usar somente o core16. HCP289 poderia ser analisado numa rede suplementar/contextual, nunca misturado silenciosamente ao core16. A etapa foi concluída conforme documentado na Seção 7; a figura-resumo agora corresponde ao painel global → Legacy28 → Core16.

## 7. Redes de interação proteína-proteína no STRING

### Objetivo e desenho analítico

A etapa STRING foi realizada após a definição das listas gênicas e teve finalidade de visualização e interpretação de relações proteína-proteína; ela não participou da definição do core16, não repriorizou genes e não constitui validação genética independente. Foram construídas três redes separadas para *Homo sapiens*, sem adição de interatores externos:

| Rede | Lista de entrada | Papel analítico |
|---|---:|---|
| `core16` | 16 genes *protein-coding* priorizados | Rede principal, derivada da regra multiestudo e multimodalidade. |
| `legacy28` | 28 genes da análise histórica Jansen–Watanabe | Sensibilidade histórica; não é independente do conjunto principal. |
| `global_gwas_catalog` | 786 entradas gênicas extraídas das associações de insônia do GWAS Catalog | Camada descritiva ampla; não é uma lista causal, nem redefine o core16. |

As redes foram configuradas como **full STRING network**, com escore de confiança alto. Na integração quantitativa, somente interações com *combined score* >= 0,700 foram contabilizadas. Os termos funcionais foram registrados quando o FDR fornecido pelo STRING era <= 0,01. Todos os canais de evidência disponíveis no STRING foram mantidos ativos, e os nós desconectados foram preservados por serem informativos sobre a ausência de um módulo proteico denso.

### Arquivos de entrada e resultados brutos

As listas submetidas e os resultados exportados permanecem em:

```text
projeto_insonia_tcc/
├── string_inputs/
│   └── core16_string_symbols.txt
├── enrichment_inputs/
│   ├── legacy_jansen_watanabe_28/
│   │   └── legacy_jansen_watanabe_28_gene_symbols.txt
│   └── gwas_catalog_insomnia_global/
│       └── gwas_catalog_insomnia_global_ensembl_ids.txt
└── Resultados_STRING/
    ├── STRING_core16/
    ├── STRING_legacy28/
    └── STRING_global_gwas_catalog/
```

O arquivo `core16_string_symbols.txt` é a lista preservada especificamente para STRING. Para Legacy28, foi usada a lista de símbolos em `enrichment_inputs/legacy_jansen_watanabe_28/`. Para a camada global, foram usados os 786 genes Ensembl do GWAS Catalog em `enrichment_inputs/gwas_catalog_insomnia_global/`, convertidos e confirmados pelo mapeamento do STRING. As tabelas `string_mapping_*.xlsx` registram os identificadores efetivamente reconhecidos pela plataforma.

Cada pasta de resultado preserva a imagem de rede, o arquivo SVG, a tabela de mapeamento de proteínas, a tabela de interações, os graus dos nós, anotações funcionais e, quando disponível, exportações de enriquecimento. Capturas de tela foram preservadas em `PRINTS_da_tela/` apenas para auditoria visual das configurações e resultados.

### Integração reprodutível e script associado

O script `scripts/25_integrate_string_results.js` consolida as três exportações manuais. Ele lê tabelas de interações e graus dos nós, deduplica arestas não direcionadas, restringe a contagem a *combined score* >= 0,700 e produz:

- `Resultados/analise_padronizada/51_string_network_summary.tsv` — resumo de entrada, mapeamento, conectividade e PPI enrichment;
- `Resultados/analise_padronizada/52_string_top_nodes.tsv` — nós ordenados por grau;
- `Resultados/analise_padronizada/53_string_key_functional_terms.tsv` — termos funcionais STRING com FDR <= 0,01;
- `Resultados/analise_padronizada/54_string_interpretation_notes.tsv` — notas de interpretação e limites;
- `Resultados/analise_padronizada/DOCUMENTACAO_STRING_E_INTERPRETACAO.md` — registro técnico resumido da etapa.

O script requer Node.js >=18 e não usa pacotes externos. Ele é uma fonte final para as tabelas integradas do STRING; não modifica as exportações originais da plataforma.

### Resultados quantitativos

| Rede | Proteínas mapeadas | Arestas >=0,700 | Nós conectados | Nós isolados | PPI enrichment |
|---|---:|---:|---:|---:|---:|
| Core16 | 16/16 | 1 | 2 | 14 | *P* = 0,0573 |
| Legacy28 | 27/28 | 2 | 4 | 23 | *P* = 0,0153 |
| GWAS Catalog global | 643/786 | 359 | 297 | 346 | Não disponível na exportação |

Na rede principal, a única interação de alta confiança foi **HEXIM1–HEXIM2** (*combined score* = 0,995). O valor de PPI enrichment de 0,0573 não sustenta conectividade global acima do acaso ao limiar convencional de 0,05; assim, a rede Core16 deve ser apresentada como conjunto heterogêneo de genes priorizados, e não como um único módulo proteico.

A Legacy28 apresentou, além de HEXIM1–HEXIM2, a interação **BTBD9–MEIS1** (*combined score* = 0,765). Seu *P* de enriquecimento PPI de 0,0153 indica mais conectividade que a esperada ao acaso, porém a rede ainda é muito esparsa e mantém caráter de sensibilidade histórica.

A rede global foi a única com escala suficiente para uma leitura de processos. Os termos mais robustos incluíram modulação da transmissão química sináptica, sinalização trans-sináptica, desenvolvimento do sistema nervoso, sinapse e compartimento somatodendrítico. Esses resultados descrevem processos compartilhados entre genes anotados em associações publicadas de insônia e devem ser interpretados como contexto biológico amplo, não como vias causalmente demonstradas para insônia.

O PPI enrichment da rede global foi registrado como “não disponível” porque essa estatística não foi incluída nas exportações preservadas da sessão STRING. Isto não significa ausência de significância e não deve ser substituído por um cálculo ad hoc; o teste depende do modelo de fundo interno do STRING.

### Regras de interpretação e apresentação no TCC

O painel de redes deve ser apresentado na ordem **GWAS Catalog global -> Legacy28 -> Core16**, visualizando o funil de evidência: associações publicadas amplas, seleção histórica e priorização multiestudo. Essa ordem não equivale a uma hierarquia de “verdade” das redes; ela mostra graus distintos de abrangência e rigor de seleção.

Uma aresta STRING representa evidência integrada de associação funcional/proteica e não prova interação física direta, regulação causal nem relação específica com insônia. Da mesma forma, anotações de fenótipo como “insônia” ou “sono” não devem ser apresentadas como descoberta funcional, pois seriam circulares em uma lista selecionada por associações de insônia.

Para o corpo do TCC, utilizar o painel das três redes e a tabela integrada `51_string_network_summary.tsv`. As tabelas de interações completas, mapeamentos e enriquecimento funcional devem permanecer como material suplementar e no repositório do pipeline.
