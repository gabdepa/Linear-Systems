#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

// GRR20197155 -> Gabriel Razzolini Pires De Paula -> grpp19

/*!
 * @brief
 *
 * @param tam Tamanho do Sistema Linear
 * @param t_egp  - Tempo da Eliminacao de Gauss com Pivoteamento Parcial(EGP)
 * @param norma_egp  - Norma da Eliminacao de Gauss com Pivoteamento Parcial(EGP) 
 * @param t_gs  Tempo do Gauss-Seidel
 * @param it_gs - Iterações do Gauss-Seidel
 * @param norma_gs Norma do Gauss-Seidel
 * @param t_ref - Tempo do Refinamento
 * @param it_ref - Iterações do Refinamento
 * @param norma_ref - Norma do Refinamento
 */
void geraTabelaResult(int tam, double t_egp, real_t norma_egp, double t_gs, int it_gs, real_t norma_gs, double t_ref, int it_ref, real_t norma_ref, int cod_return)
{
  printf("n:\t| t_egp:\t| normaResiduo_egp:\t| t_gs:\t\t\t| it_gs:\t| normaResiduo_gs:\t| t_ref:\t\t| it_ref:\t| normaResiduo_ref:\n");
  printf("%d\t| %g\t| %10g\t\t| %g\t\t| %d\t\t| %10g\t\t| %g\t\t| %d\t\t| %10g\n",
         tam,
         t_egp,
         norma_egp,
         t_gs,
         it_gs,
         norma_gs,
         t_ref,
         it_ref,
         norma_ref);
}

int main()
{
  // inicializa gerador de números aleatóreos
  srand(202202);

  // Sistema Linear
  SistLinear_t *SL;

  // Código de Retorno das funções
  int cod_return;

  // Tamanho das Matrizes
  int tamanhos[10] = {10, 30, 50, 128, 256, 512, 1000, 2000, 3000};

  // Laço que percorre os 3 tipos de Matriz: Genérica, Hilbert e Diagonalmente Dominate
  for (int tipoMat = 0; tipoMat < 3; ++tipoMat)
  {
    printf("\n------- Tipo Matriz %d -------\n", tipoMat);

    // Laço que percorre os tamanhos para cada Matriz do vetor tamanhos[10]
    for (int i = 0; i < 9; ++i)
    {
      // Atualiza variável do tamanho para receber valor de tamanhos[i]
      int tamAtual = tamanhos[i];

      // Inicializa variáveis com vetor de soluções para EGP, Gauss-Seidel e Refinamento, respectivamente
      real_t sol_egp[tamAtual + 1], sol_gs[tamAtual + 1], sol_ref[tamAtual + 1];

      // Inicializa variável que guarda valor do número de iterações para Refinamento e Gauss-Seidel
      int it_ref, it_gs;
      it_ref = it_gs = 0;

      // Inicializa variável que guarda tempo gasto nas funções EGP, Gauss-Seidel e Refinamento, respectivamente
      double t_egp, t_gs, t_ref;
      t_egp = t_gs = t_ref = 0;

      // Alocação do sistema Linear com tamanho atual
      SL = alocaSisLin(tamAtual);

      // Inicialização do sistema Linear com tipo da Matriz e coeficiente máximo que cada coeficiente pode assumir
      iniSisLin(SL, tipoMat, COEF_MAX);

      // Laço que zera todos os valores do vetor de soluções de EGP, Gauss-Seidel e Refinamento
      for (int i = 0; i < tamAtual; ++i)
      {
        sol_egp[i] = sol_gs[i] = sol_ref[i] = 0;
      }
      
      // Guarda código de retorno da função Eliminação de Gauss com Pivoteamento Parcial(EGP)
      cod_return = eliminacaoGauss(SL, sol_egp, &t_egp);

      // Guarda valor do número de iterações feitas pela função Gauss-Seidel
      it_gs = gaussSeidel(SL, sol_gs, ERRO, &t_gs);

      // Inicializa vetor solucao do refinamento com a solucao do EGP
      for (int i = 0; i < tamAtual; ++i)
      {
        sol_ref[i] = sol_gs[i];
      }

      // Guarda valor do número de iterações feitas pela função de Refinamento
      it_ref = refinamento(SL, sol_ref, ERRO, &t_ref);

      // Gera tabela na stdout com os resultados de cada função
      geraTabelaResult(tamAtual, t_egp, normaL2Residuo(SL, sol_egp), t_gs, it_gs, normaL2Residuo(SL, sol_gs), t_ref, it_ref, normaL2Residuo(SL, sol_gs), cod_return);
      
      // Liberação de memória do Sistema Alocado
      liberaSisLin(SL);
    }
  }

  return 0;
}
