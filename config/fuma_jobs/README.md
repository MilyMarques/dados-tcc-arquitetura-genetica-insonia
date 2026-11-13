# Parâmetros e metadados dos jobs FUMA

`snp2gene/` contém o arquivo `params.config` exportado de cada job SNP2GENE usado no projeto: quatro estudos principais, dois proxies OpenGWAS e FinnGen R13 corrigido.

`gene2func/` contém os parâmetros exportados dos jobs GENE2FUNC das listas Core16, Legacy28, HCP289, FinnGen R13 e GWAS Catalog global.

`submission_metadata/` preserva o build genômico, esquema de colunas, alelo de efeito e outras informações de preparação dos inputs submetidos aos jobs SNP2GENE.

Os ZIPs completos do FUMA e as estatísticas-resumo não são distribuídos neste repositório. Consulte `manifest.json` para a relação de arquivos preservados.

Os arquivos preservados correspondem aos jobs finais que alimentam as tabelas integradas do release. Configurações históricas de Jansen e Watanabe, anteriores à reexecução final, não são utilizadas como fonte dos resultados deste repositório.
