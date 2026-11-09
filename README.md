# Arquitetura genética da insônia: pipeline reprodutível

Este repositório contém a versão organizada do pipeline do TCC: preparação/QC de sumstats, FUMA SNP2GENE e MAGMA, integração multiestudo, enriquecimento funcional e redes STRING.

## Resultado principal

O conjunto **Core16** reúne genes *protein-coding* com significância MAGMA por Bonferroni em pelo menos duas **unidades analíticas** e ao menos uma modalidade de mapeamento SNP2GENE. Ele está em `results/tables/22_core16_final_evidence_table.tsv`. Jansen e Watanabe foram tratados como uma única unidade UK Biobank na regra de consenso devido à sobreposição amostral; portanto, a regra representa convergência entre análises e não independência completa de participantes.

## Organização

- `scripts_em_python/`: conversão, QC e validação de inputs FUMA.
- `scripts_em_python_e_node/`: integrações complementares, incluindo GWAS Catalog, FinnGen corrigido e STRING.
- `scripts_em_R/`: integração dos resultados FUMA/MAGMA.
- `inputs/`: listas pequenas reutilizáveis para análises de sensibilidade e STRING.
- `results/tables/`: tabelas finais para o manuscrito.
- `results/figures/`: figuras STRING selecionadas.
- `results/string_exports/`: exportações tabulares originais do STRING que sustentam a integração.
- `docs/`: metodologia, interpretação e guia de escrita.
- `data_availability/`: política para dados brutos não versionados.

## Como interpretar as redes

As três redes correspondem a três níveis de evidência: Core16 (principal), Legacy28 (sensibilidade histórica) e GWAS Catalog global (camada descritiva). A última revela convergência para sinalização sináptica, transmissão trans-sináptica, desenvolvimento do sistema nervoso e adesão celular, mas não redefine o Core16 nem demonstra causalidade. Veja `results/tables/51_string_network_summary.tsv` e `docs/DOCUMENTACAO_STRING_E_INTERPRETACAO.md`.

## Dados não incluídos

Não suba estatísticas-resumo brutas, dados individuais, arquivos VCF/BCF, bibliotecas R, ZIPs completos do FUMA ou screenshots de interface. Eles permanecem no diretório de trabalho local e devem ser recuperados nas fontes oficiais, conforme licenças e versões descritas em `docs/`.

## Reprodutibilidade

Esta versão é um **release de código, parâmetros e resultados derivados selecionados**. Para reproduzir o fluxo integral, obtenha as estatísticas-resumo e os ZIPs dos jobs FUMA nas fontes originais, respeitando as licenças, e reconstrua a estrutura de diretórios descrita em `docs/DOCUMENTACAO_COMPLETA_METODOLOGIA_E_SCRIPTS.md`. O FinnGen final é exclusivamente o job corrigido 763972. Resultados FinnGen anteriores são históricos e não devem ser usados como evidência de replicação.

## Antes de tornar público

Revise autoria, afiliações, referências, licença escolhida pelo grupo e se as políticas de cada repositório permitem redistribuir cada arquivo derivado.
