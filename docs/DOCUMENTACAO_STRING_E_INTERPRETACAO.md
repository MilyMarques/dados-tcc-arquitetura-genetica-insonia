# STRING: integração e interpretação

## Configuração comum

As três redes foram geradas manualmente no STRING para *Homo sapiens*, tipo **full STRING network**, sem nós adicionais, com escore de confiança alto (>=0,700) e controle de FDR estrito (1%). O arquivo `51_string_network_summary.tsv` é a síntese quantitativa; `52_string_top_nodes.tsv` lista os nós de maior grau; `53_string_key_functional_terms.tsv` preserva os termos com FDR <=0,01.

## Leitura por camada

- **Core16:** 16 proteínas mapeadas, 1 aresta (HEXIM1–HEXIM2; score 0,995), PPI enrichment *p*=0,0573. A figura mostra heterogeneidade funcional entre os genes priorizados; não sustenta um único módulo proteico. É a figura de rede principal porque os genes foram definidos com a regra mais rigorosa.
- **Legacy28:** 27/28 genes mapeados (o transcrito não codificante RP11-463D19.2 não possui proteína humana correspondente), 2 arestas e PPI enrichment *p*=0,0153. A conexão adicional é BTBD9–MEIS1 (score 0,765). Trate como análise de sensibilidade histórica, não como validação independente.
- **Global GWAS Catalog:** 643/786 entradas foram mapeadas a proteínas, com 359 arestas. Esta escala revela a arquitetura funcional descritiva de genes associados à insônia em estudos publicados: sinalização sináptica, transmissão trans-sináptica, desenvolvimento do sistema nervoso e adesão celular. A imagem é visualmente informativa, mas não pode redefinir o Core16.

## Como reportar

Use o painel das três redes para explicar o funil analítico: ampla evidência publicada (global) -> seleção histórica (Legacy28) -> conjunto replicado e priorizado (Core16). Não conclua que a rede global prova mecanismo causal ou que o PPI da Legacy28 represente uma via única.
