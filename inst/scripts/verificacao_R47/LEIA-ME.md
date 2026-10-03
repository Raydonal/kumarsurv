# verificacao_R47: scripts de auditoria da Rodada 47

Estes scripts sustentam os achados da auditoria da Rodada 47 (relatório `RELATORIO_RODADA47.md`,
que acompanha o manuscrito). Eles
**não** geram tabelas do artigo: as tabelas vêm dos scripts registrados em `scripts/run_all.R`.
Estão aqui para que qualquer coautor possa refazer as verificações.

Para rodar todos, na ordem das dependências, a partir da raiz do pacote (cerca de 15 min):

```bash
LC_ALL=C.utf8 Rscript scripts/verificacao_R47/rodar_verificacao_R47.R
```

As transcrições vão para `resultados/verificacao_R47/*.txt`. Os artefatos `R47_*.rds` vão para
`resultados/`.

| script | o que verifica | achado |
|---|---|---|
| `R47_teste_kalman.R` | forma quadrática (Kalman) e log-determinante contra os valores exatos (Durbin–Levinson completo), na aplicação e no desenho de Varin | longe da fronteira, as vias coincidem em todos os dígitos impressos |
| `R47_teste_hessiana.R` | EPs com o passo padrão do `numDeriv` (`d = 0.1`), com `d = 1e-3` e com `optimHess` | o passo padrão leva o AR para fora da região estacionária: EP(ψ₁) = NaN na aplicação e EPs de dependência subestimados no desenho de Varin (item 2.1) |
| `R47_otimo_global.R` | se o ajuste adotado (Kumaraswamy, ARMA(1,1)) é o máximo global | o submodelo aninhado com forma constante tem log-verossimilhança maior (2543.21 contra 2542.27) |
| `R47_perfil_ar1_global.R` | multipartida e perfil de ψ₁ até 0.999, nas duas bacias | a verossimilhança cresce em direção à fronteira (cerca de 2545.4 em 0.999) |
| `R47_vero_exata.R` | `nll_copula` contra a verossimilhança exata junto à fronteira; `kal()` reconstrói o cálculo anterior (2000 pesos ψ) | com a biblioteca corrigida, coincidem em todos os pontos, inclusive junto à fronteira |
| `R47_global_exata.R` | máximo e perfil com a verossimilhança exata, restrita a \|AR\| ≤ 0.997 | a verossimilhança da Kumaraswamy cresce até a fronteira; a da beta tem máximo interior |
| `R47_valida_lib.R` | biblioteca corrigida, anterior (`copula_ml_pre_R47.R`) e exata, variando só o AR, e no ponto de maior verossimilhança que a biblioteca anterior encontrava (ψ₁ = 0.9995) | a corrigida coincide com a exata em todos os pontos; a anterior diverge a partir de 0.996 e, naquele ponto, dá 2562.99 contra 2533.04 da exata |
| `R47_periodo_modos.R` | em que modo (AR interior ou de fronteira) caiu cada ajuste da comparação de periodicidade do v11 | o semestral caiu no modo de fronteira (0.994) e os demais no interior. A coluna "melhor" usa penalidade fora da faixa e subestima o máximo restrito; os valores da Tabela S12 vêm de `perfil_ar1_exato.R` |
| `R47_faseI.R` | cobertura dos limites na Fase I | base da coluna Fase I da Tabela 6 |
| `R47_calib_controle.R` | primeira versão da calibração sob controle, com os limites marginais e condicionais | origem do `calib_controle_SEM.R` registrado |
| `R47_calib_robustez.R` | calibração sob controle excluindo só as semanas do Farrington, ou também as da CUSUM | a conclusão não depende da definição do surto |
| `R47_ep_varin.R` | EPs do desenho de Varin recalculados com `d = 1e-3` nas estimativas **arquivadas**, regenerando os dados de cada réplica pela mesma semente | `VERIFICA=1` confere que o reajuste reproduz as estimativas arquivadas |
| `R47_tab_varin.R` | Quadro 8 a partir dos EPs acima | a razão DPr/EP fica entre 0.91 e 1.13; três réplicas atípicas (ver abaixo) |

`copula_ml_pre_R47.R` é a biblioteca como estava antes da Rodada 47 (a mesma de
`scripts_BACKUP_pre_rodada47/copula_ml.R`). Só `R47_valida_lib.R` a usa, para a comparação.

Com a biblioteca corrigida, `14_mc_validacao_B300.R` grava os EPs por réplica e a coluna DPr.
O Quadro 8 do suplementar sai dele, e `R47_ep_varin.R`/`R47_tab_varin.R` passam a ser uma
conferência independente: avaliadas nas mesmas estimativas, as duas vias dão EPs idênticos, bit a
bit. Durante a auditoria, `R47_ep_varin.R` foi avaliado nas estimativas arquivadas antes da
correção, que diferem das atuais por até 0.0098 em 23 réplicas e por menos de 1e-6 nas demais.
Com isso, a média dos EPs mudava na terceira casa (0.162 contra 0.160 no intercepto da forma),
porque em poucas réplicas a hessiana é mal condicionada. A razão DPr/EP ia de 0.91 a 1.11
naquela avaliação e vai de 0.91 a 1.13 na atual.
