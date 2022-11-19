#ifndef __METODOS_H__
#define __METODOS_H__

// Parâmetros para teste de convergência
#define MAXIT 50     // número máximo de iterações para métodos iterativos
#define ERRO 1.0e-10  // Tolerância para critérios de parada em métodos iterativos

/*!
  \brief Encontra maior valor do coeficiente da linha da Matriz

  \param SL Ponteiro para o sistema linear
  \param i Linha indicada para procura

  \return Índice do maior coeficiente
*/ 
int encontraMax(SistLinear_t *SL, int i);


/*!
  \brief Troca linha do maior coeficiente com linha indicada por i

  \param SL Ponteiro para o sistema linear
  \param i Linha indicada para troca
  \param max_index Índice do maior coeficiente A

*/
void trocaLinha(SistLinear_t *SL, int i, int max_index);


/*!
  \brief Método da Eliminação de Gauss com Pivoteamento Parcial(EGP)

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. 0 em caso de sucesso.
*/int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal);


/*!
  \brief Calcula Norma Vetorial do Erro em Absoluto

  \param aux Ponteiro para o vetor solução calculado anteriormente
  \param x Ponteiro para o vetor solução atual
  \param tam tamanho do vetor de soluções

  \return Norma calculada
*/
real_t calcularNormaMaxErroAbsoluto(real_t *aux, real_t *x, int tam);


/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
*/
int gaussSeidel (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal);



/*!
  \brief Calcula a multiplicação de A*X da matriz

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param solution Solução

*/void multiMatrix(SistLinear_t *SL, real_t *x, real_t *solution);


/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear

  \return Norma L2 do resíduo
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *x);

/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
int refinamento (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal);

#endif // __METODOS_H__

