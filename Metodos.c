#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"


int encontraMax(SistLinear_t *SL, int i)
{
  // guarda o maior valor como o primeiro elemento da linha
  real_t max = SL->A[i][i];

  // guarda o índice do maior
  int max_index = i;

  // loop para percorrer as linhas da matriz
  for (int l = i; l < SL->n; ++l)
  {
    // Verifica se o elemento encontrado é maior que o maior até então
    if (ABS(SL->A[l][i]) > max)
    {
      max = SL->A[l][i];
      max_index = l;
    }
  }
  return max_index;
}

void trocaLinha(SistLinear_t *SL, int i, int max_index)
{
  // Percorre todos os elementos da linha, mudando as colunas 'k'
  for (int k = 0; k <= SL->n; ++k)
  {
    real_t temp = SL->A[i][k];
    SL->A[i][k] = SL->A[max_index][k];
    SL->A[max_index][k] = temp;
  }

  // Troca o termo independente do vetor de termos independentes
  real_t temp = SL->b[max_index];
  SL->b[max_index] = SL->b[i];
  SL->b[i] = temp;
}

int eliminacaoGauss(SistLinear_t *SL, real_t *x, double *tTotal)
{
  // inicializa tempo da função
  double time = timestamp();

  /* para cada linha a partir da primeira */
  for (int i = 0; i < SL->n; ++i)
  {
    unsigned int max_index = encontraMax(SL, i);
    if (i != max_index)
      trocaLinha(SL, i, max_index);

    for (int k = i + 1; k < SL->n; ++k)
    {
      real_t m = SL->A[k][i] / SL->A[i][i]; // pivô

      // Verificação se resultou em NaN ou +/- infinito
      if(isnan(m) || isinf(m))
      {
        fprintf(stderr, "Erro gauss - linha %i m: %g é NaN ou +/-Infinito\n", i, m);
        return -1;
      }
      // zera o elemento para diminuir o número de operações
      SL->A[k][i] = 0.0;

      for (int j = i + 1; j < SL->n; ++j)
        SL->A[k][j] -= SL->A[i][j] * m; //
      x[k] -= x[i] * m;
    }

    // retrosubstituição
    for (int i = SL->n - 1; i >= 0; --i)
    {
      real_t sum = 0;
      for (int j = i; j < SL->n; ++j)
      {
        sum = sum + SL->A[i][j] * SL->b[j];
      }
      x[i] = (SL->A[i][SL->n + 1] - sum) / SL->A[i][i];
    }
  }

  // Coloca tempo decorrido para execução da função no vetor tTotal do programa e libera matriz anteriormente alocada;
  time = timestamp() - time;
  tTotal[0] = time;

  return 1;
}

real_t calcularNormaMaxErroAbsoluto(real_t *aux, real_t *x, int tam)
{
  real_t max = ABS(x[0] - aux[0]);

  // Percorre todos os termos
  for (int i = 1; i < tam; ++i)
  {
    // Se o módulo da subtração é maior que max
    if (ABS(x[i] - aux[i]) > max)
    {
      max = ABS(x[i] - aux[i]);
    }
  }
  return max;
}

int gaussSeidel(SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{
  tTotal[0] = timestamp();

  int num_iterations = 0;

  real_t error_calculated;
  real_t *x_calculated;
  x_calculated = (real_t *)calloc(SL->n, sizeof(real_t));

  // For que preenche o vetor solução com o chute inicial, "0.00"
  for (int i = 0; i < SL->n; ++i)
  {
    x[i] = 0.00;
  }

  // gauss-Seidel method
  while (1)
  {
    for (int i = 0; i < SL->n; ++i)
    {
      for (int j = 0; j < SL->n; j++)
      {
        if (j < i)
        {
          // substituição do x calculado para utilização na próxima linha
          x_calculated[i] += -(SL->A[i][j]) * x_calculated[j];
        }
        else if (j > i)
        {
          // linha que calcula a soma dos coeficientes multiplicado pelo chute inicial
          x_calculated[i] += -(SL->A[i][j]) * x[j];
        }
      }
      // linha que calcula o valor de c(n) - sum / a(n,n)
      x_calculated[i] += SL->b[i];
      x_calculated[i] /= SL->A[i][i];

      // Verificação se resultou em NaN ou +/- infinito
      if(isnan(x_calculated[i]) || isinf(x_calculated[i]))
      {
        fprintf(stderr, "Erro gauss-Seidel - linha %i x_calculated[i]: %g é NaN ou +/-Infinito\n", i, x_calculated[i]);
        return -1;
      }
    }
    // Incrementa número de iterações
    ++num_iterations;

    // calcula erro em relação ao vetor de solução
    error_calculated = calcularNormaMaxErroAbsoluto(x_calculated, x, SL->n);
    // Critério de Parada
    if ((error_calculated <= erro) || (num_iterations >= MAXIT))
    {
      for (int i = 0; i < SL->n; ++i)
        x[i] = x_calculated[i];

      tTotal[0] = timestamp() - tTotal[0];

      // Retorna número de iterações
      return num_iterations;
    }

    // else...
    for (int i = 0; i < SL->n; ++i)
    {
      x[i] = x_calculated[i];
      x_calculated[i] = (real_t)0;
    }
  }
}

void multiMatrix(SistLinear_t *SL, real_t *x, real_t *solution)
{
  // solution = [A]*[X]
  for (int i = 0; i < SL->n; ++i)
    solution[i] = 0;

  for (int i = 0; i < SL->n; ++i)
  {
    for (int j = 0; j < SL->n; ++j)
      solution[i] += SL->A[i][j] * x[j];
  }
}

real_t normaL2Residuo(SistLinear_t *SL, real_t *x)
{
  real_t norma = 0;
  real_t residue[SL->n+1]; 
  // A*X
  real_t matrix_solution[SL->n+1];
  // alocação de memória para variável Auxiliar que guarda o valor da multiplicação da matriz
  // matrix_solution = (real_t *)calloc(SL->n, sizeof(real_t));

  // função que calcula a Multiplicação A*X da matriz
  multiMatrix(SL, x, matrix_solution);

  // r0 = b0 - x0, r1 = b1 - x1, r2 = b2 - x2....
  for (int i = 0; i < SL->n; i++)
    residue[i] = SL->b[i] - matrix_solution[i];

  // Calcula da norma -> norma = sqrt(r1²+r2²+r3²+...)
  for (int i = 0; i < SL->n; i++)
    norma += residue[i] * residue[i];
  norma = sqrt(norma);

  // Verificação se resultou em NaN ou +/- infinito
  if(isnan(norma) || isinf(norma))
  {
    fprintf(stderr, "Erro norma: %g é NaN ou +/-Infinito\n", norma);
    return -1;
  }

  return norma;
}

int refinamento(SistLinear_t *SL, real_t *x, real_t erro, double *tTotal)
{
  // inicializa tempo da função
  double time = timestamp();

  // inicializa contador de iterações
  int iteracao = 0;

  // inicializa variáveis necessárias para que seja realizado refinamento
  real_t *aux, *w, *r, norma;

  // alocação de memória para variável auxiliar "aux"
  aux = (real_t *)calloc(SL->n, sizeof(real_t));


  // alocação de memória para variável "w" erro
  w = (real_t *)calloc(SL->n, sizeof(real_t));

  // alocação de memória para variável "resíduo" atual
  r = (real_t *)calloc(SL->n, sizeof(real_t));

  // guarda valor da Norma calculado pela função
  norma = normaL2Residuo(SL, x);

  // copia sistema linear original para nova iteração do refinamento
  SistLinear_t *newSL = alocaSisLin(SL->n);

  while ((norma > 5.0) || (iteracao < MAXIT))
  {
    iteracao++;

    multiMatrix(SL, x, aux);

    for (int i = 0; i < SL->n; ++i)
      r[i] = SL->b[i] - aux[i];

    // guarda SL antigo em nova variável para atualização do resíduo
    for (int i = 0; i < SL->n; ++i)
      for (int j = 0; j < SL->n; ++j)
        newSL->A[i][j] = SL->A[i][j];

    // vetor de termos indepentes que recebe valor do resíduo[i]
    for (int i = 0; i < SL->n; ++i)
      newSL->b[i] = r[i];

    // ajusta variável tamanho
    newSL->n = SL->n;

    // Chamada da eliminação de Gauss que resolve-> A*w = r
    eliminacaoGauss(newSL, w, tTotal);
    // copia a nova solucao obtida para novo u(u_ant+1)
    for (int i = 0; i < SL->n; ++i)
    {
      x[i] += w[i];
      w[i] = 0;
    }

    norma = normaL2Residuo(SL, x);

    if(norma == -1)
    {
      time = timestamp() - time;
      tTotal[0] = time;
      return -1;
    } 
  }

  // Liberação de estruturas utilizadas
  liberaSisLin(newSL);
  free(aux);
  free(w);
  free(r);

  // Coloca tempo decorrido para execução da função no vetor tTotal do programa
  time = timestamp() - time;
  tTotal[0] = time;

  return iteracao;
}