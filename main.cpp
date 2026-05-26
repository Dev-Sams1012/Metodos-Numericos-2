// =============================================================================
// main.cpp
// Interface de linha de comando (CLI) para explorar todos os métodos numéricos.
//
// Estrutura de menus:
//   [1] Diferenciação Numérica
//       [1] Derivada Primeira (Forward/Backward/Central)
//       [2] Derivada Segunda (Central)
//       [3] Derivada Terceira (Forward/Central)
//   [2] Integração Numérica
//       [1] Newton-Cotes Fechado (Trapézio, Simpson 1/3, Simpson 3/8)
//       [2] Newton-Cotes Aberto (Trapézio Aberto, Milne)
//       [3] Gauss-Legendre (n = 2..5)
//       [4] Gauss Especial (Hermite, Laguerre, Chebyshev)
//       [5] Singularidades (Exponencial Simples / Dupla)
//   [3] Autovalores e Autovetores
//       [1] Método da Potência (Regular / Inverso / Deslocamento)
//       [2] Transformação de Householder (tridiagonalização)
//       [3] Método de Jacobi
//       [4] Método QR
//       [5] SVD
//   [0] Sair
// =============================================================================

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <windows.h>

// Core
#include "core/FunctionParser.hpp"

// Diferenciação
#include "differentiation/FirstDerivative.hpp"
#include "differentiation/SecondDerivative.hpp"
#include "differentiation/ThirdDerivative.hpp"

// Integração
#include "integration/NewtonCotes.hpp"
#include "integration/GaussLegendre.hpp"
#include "integration/SpecialGauss.hpp"
#include "integration/SingularityStrategies.hpp"

// Autovalores
#include "eigenvalues/PowerMethod.hpp"
#include "eigenvalues/HouseholderMethod.hpp"
#include "eigenvalues/JacobiMethod.hpp"
#include "eigenvalues/QRMethod.hpp"
#include "eigenvalues/SVDMethod.hpp"

// Utils
#include "Matrix.hpp"

using namespace nm;

// =============================================================================
// Helpers de UI
// =============================================================================

static void clearLine() {
    std::cout << std::string(70, '-') << "\n";
}

static void header(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(70, '=') << "\n";
}

static void subheader(const std::string& title) {
    std::cout << "\n" << std::string(50, '-') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(50, '-') << "\n";
}

/// Lê uma linha e converte para T; repete se inválido.
template<typename T>
T readValue(const std::string& prompt) {
    T value{};
    while (true) {
        std::cout << prompt;
        std::string line;
        std::getline(std::cin, line);
        std::istringstream iss(line);
        if (iss >> value) break;
        std::cout << "  [!] Entrada inválida. Tente novamente.\n";
    }
    return value;
}

/// Lê inteiro dentro de [lo, hi].
static int readInt(const std::string& prompt, int lo, int hi) {
    while (true) {
        int v = readValue<int>(prompt);
        if (v >= lo && v <= hi) return v;
        std::cout << "  [!] Valor fora do intervalo [" << lo << ", " << hi << "].\n";
    }
}

/// Lê expressão e monta FunctionParser 1D.
static std::shared_ptr<FunctionParser> readFunction1D() {
    auto fp = std::make_shared<FunctionParser>();
    while (true) {
        std::cout << "  f(x) = ";
        std::string expr;
        std::getline(std::cin, expr);
        try {
            fp->parse(expr);
            std::cout << "  → Função aceita: f(x) = " << expr << "\n";
            return fp;
        } catch (const std::exception& e) {
            std::cout << "  [!] " << e.what() << "\n";
        }
    }
}

/// Lê expressão e monta FunctionParser 2D.
static std::shared_ptr<FunctionParser> readFunction2D() {
    auto fp = std::make_shared<FunctionParser>();
    while (true) {
        std::cout << "  f(x,y) = ";
        std::string expr;
        std::getline(std::cin, expr);
        try {
            fp->parse2D(expr);
            std::cout << "  → Função aceita: f(x,y) = " << expr << "\n";
            return fp;
        } catch (const std::exception& e) {
            std::cout << "  [!] " << e.what() << "\n";
        }
    }
}

/// Lê uma matriz n×m do usuário linha por linha.
static Matrix readMatrix(int rows, int cols) {
    Matrix M(rows, cols);
    std::cout << "  Digite os elementos linha a linha (separados por espaço):\n";
    for (int i = 0; i < rows; ++i) {
        while (true) {
            std::cout << "  Linha " << i + 1 << ": ";
            std::string line;
            std::getline(std::cin, line);
            std::istringstream iss(line);
            bool ok = true;
            for (int j = 0; j < cols; ++j) {
                double v;
                if (!(iss >> v)) { ok = false; break; }
                M(i, j) = v;
            }
            if (ok) break;
            std::cout << "  [!] Esperados " << cols << " valores. Tente novamente.\n";
        }
    }
    return M;
}

/// Imprime vetor.
static void printVec(const std::string& label, const std::vector<double>& v, int prec = 8) {
    std::cout << "  " << label << " = [";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << std::fixed << std::setprecision(prec) << v[i];
        if (i + 1 < v.size()) std::cout << ",  ";
    }
    std::cout << "]\n";
}

/// Imprime resultado escalar.
static void printResult(const std::string& method, const std::string& label,
                        double value, int prec = 10) {
    clearLine();
    std::cout << "  Método  : " << method << "\n";
    std::cout << "  " << label << " = "
              << std::fixed << std::setprecision(prec) << value << "\n";
    clearLine();
}

// =============================================================================
// MÓDULO 1 — DIFERENCIAÇÃO NUMÉRICA
// =============================================================================

static void menuPrimeiraDerivada() {
    subheader("Derivada Primeira f'(x)");
    std::cout << "  Métodos disponíveis:\n";
    std::cout << "  [1] Forward  2 pontos  O(h)\n";
    std::cout << "  [2] Forward  3 pontos  O(h²)\n";
    std::cout << "  [3] Forward  4 pontos  O(h³)\n";
    std::cout << "  [4] Backward 2 pontos  O(h)\n";
    std::cout << "  [5] Backward 3 pontos  O(h²)\n";
    std::cout << "  [6] Backward 4 pontos  O(h³)\n";
    std::cout << "  [7] Central  2 pontos  O(h²)\n";
    std::cout << "  [8] Central  4 pontos  O(h⁴)\n";

    int choice = readInt("\n  Escolha o método: ", 1, 8);
    FirstDerivMethod methods[] = {
        FirstDerivMethod::FORWARD_2PT,  FirstDerivMethod::FORWARD_3PT,
        FirstDerivMethod::FORWARD_4PT,  FirstDerivMethod::BACKWARD_2PT,
        FirstDerivMethod::BACKWARD_3PT, FirstDerivMethod::BACKWARD_4PT,
        FirstDerivMethod::CENTRAL_2PT,  FirstDerivMethod::CENTRAL_4PT
    };
    FirstDerivMethod m = methods[choice - 1];

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();

    double x = readValue<double>("\n  Ponto x = ");
    double h = readValue<double>("  Passo  h = ");

    FirstDerivative deriv(fp, m);
    double result = deriv.compute(x, h);
    printResult(deriv.name(), "f'(" + std::to_string(x) + ")", result);
}

static void menuSegundaDerivada() {
    subheader("Derivada Segunda f''(x)");
    std::cout << "  [1] Central 3 pontos  O(h²)\n";
    std::cout << "  [2] Central           O(h⁴)\n";

    int choice = readInt("\n  Escolha o método: ", 1, 2);
    SecondDerivMethod m = (choice == 1) ? SecondDerivMethod::CENTRAL_3PT
                                        : SecondDerivMethod::CENTRAL_4PT;

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double x = readValue<double>("\n  Ponto x = ");
    double h = readValue<double>("  Passo  h = ");

    SecondDerivative deriv(fp, m);
    double result = deriv.compute(x, h);
    printResult(deriv.name(), "f''(" + std::to_string(x) + ")", result);
}

static void menuTerceiraDerivada() {
    subheader("Derivada Terceira f'''(x)");
    std::cout << "  [1] Forward  O(h²)\n";
    std::cout << "  [2] Central  O(h⁴)\n";

    int choice = readInt("\n  Escolha o método: ", 1, 2);
    ThirdDerivMethod m = (choice == 1) ? ThirdDerivMethod::FORWARD_O2
                                       : ThirdDerivMethod::CENTRAL_O4;

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double x = readValue<double>("\n  Ponto x = ");
    double h = readValue<double>("  Passo  h = ");

    ThirdDerivative deriv(fp, m);
    double result = deriv.compute(x, h);
    printResult(deriv.name(), "f'''(" + std::to_string(x) + ")", result);
}

static void menuDiferenciacao() {
    while (true) {
        header("DIFERENCIAÇÃO NUMÉRICA");
        std::cout << "  [1] Derivada Primeira\n";
        std::cout << "  [2] Derivada Segunda\n";
        std::cout << "  [3] Derivada Terceira\n";
        std::cout << "  [0] Voltar\n";

        int op = readInt("\n  Opção: ", 0, 4);
        if (op == 0) break;
        try {
            switch (op) {
                case 1: menuPrimeiraDerivada(); break;
                case 2: menuSegundaDerivada();  break;
                case 3: menuTerceiraDerivada(); break;
            }
        } catch (const std::exception& e) {
            std::cout << "\n  [ERRO] " << e.what() << "\n";
        }
        std::cout << "\n  Pressione ENTER para continuar...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// =============================================================================
// MÓDULO 2 — INTEGRAÇÃO NUMÉRICA
// =============================================================================

static void menuNewtonCotesFechado() {
    subheader("Newton-Cotes Fechado");
    std::cout << "  [1] Regra do Trapézio\n";
    std::cout << "  [2] Simpson 1/3\n";
    std::cout << "  [3] Simpson 3/8\n";

    int choice = readInt("\n  Escolha: ", 1, 3);
    NewtonCotesMethod methods[] = {
        NewtonCotesMethod::TRAPEZOID,
        NewtonCotesMethod::SIMPSON_1_3,
        NewtonCotesMethod::SIMPSON_3_8
    };

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double a = readValue<double>("\n  Limite inferior a = ");
    double b = readValue<double>("  Limite superior b = ");

    NewtonCotes integrator(fp, methods[choice - 1]);
    double result = integrator.compute(a, b);
    printResult(integrator.name(),
                "∫[" + std::to_string(a) + "," + std::to_string(b) + "] f(x)dx",
                result);
}

static void menuNewtonCotesAberto() {
    subheader("Newton-Cotes Aberto");
    std::cout << "  [1] Trapézio Aberto\n";
    std::cout << "  [2] Regra de Milne\n";

    int choice = readInt("\n  Escolha: ", 1, 2);
    NewtonCotesMethod m = (choice == 1) ? NewtonCotesMethod::OPEN_TRAPEZOID
                                        : NewtonCotesMethod::MILNE;

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double a = readValue<double>("\n  Limite inferior a = ");
    double b = readValue<double>("  Limite superior b = ");

    NewtonCotes integrator(fp, m);
    double result = integrator.compute(a, b);
    printResult(integrator.name(),
                "∫[" + std::to_string(a) + "," + std::to_string(b) + "] f(x)dx",
                result);
}

static void menuGaussLegendre() {
    subheader("Quadratura de Gauss-Legendre");
    std::cout << "  Número de pontos n = 2, 3, 4 ou 5\n";
    int n = readInt("\n  n = ", 2, 5);

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double a = readValue<double>("\n  Limite inferior a = ");
    double b = readValue<double>("  Limite superior b = ");

    GaussLegendre integrator(fp, n);
    double result = integrator.compute(a, b);
    printResult(integrator.name(),
                "∫[" + std::to_string(a) + "," + std::to_string(b) + "] f(x)dx",
                result);
}

static void menuGaussEspecial() {
    subheader("Quadraturas Especiais de Gauss (5 pontos)");
    std::cout << "  [1] Gauss-Hermite   — ∫_{-∞}^{+∞} e^{-x²} f(x) dx\n";
    std::cout << "  [2] Gauss-Laguerre  — ∫_{0}^{+∞}  e^{-x}  f(x) dx\n";
    std::cout << "  [3] Gauss-Chebyshev — ∫_{-1}^{1}  f(x)/√(1-x²) dx\n";

    int choice = readInt("\n  Escolha: ", 1, 3);
    SpecialGaussType types[] = {
        SpecialGaussType::HERMITE,
        SpecialGaussType::LAGUERRE,
        SpecialGaussType::CHEBYSHEV
    };

    int n = 5;
    if (choice == 3) {
        std::cout << "  Número de pontos de Chebyshev (1-20): ";
        n = readInt("", 1, 20);
    }

    std::cout << "\n  ATENÇÃO: digite apenas a parte f(x) (sem o peso).\n";
    std::cout << "  Ex.: para integrar e^{-x²}·cos(x), use Hermite e escreva cos(x)\n";
    std::cout << "\n  Digite f(x):\n";
    auto fp = readFunction1D();

    SpecialGauss integrator(fp, types[choice - 1], n);
    double result = integrator.compute();
    printResult(integrator.name(), "I", result);
}

static void menuSingularidades() {
    subheader("Estratégias para Singularidades");
    std::cout << "  [1] Exponencial Simples  [x(s) = m + r·tanh(s)]\n";
    std::cout << "  [2] Exponencial Dupla    [x(s) = m + r·tanh(π/2·sinh(s))]\n";

    int choice = readInt("\n  Escolha: ", 1, 2);
    SingularityMethod m = (choice == 1) ? SingularityMethod::SINGLE_EXP
                                        : SingularityMethod::DOUBLE_EXP;

    std::cout << "\n  Digite a função f(x) (pode ter singularidade nos extremos):\n";
    auto fp = readFunction1D();
    double a = readValue<double>("\n  Limite inferior a = ");
    double b = readValue<double>("  Limite superior b = ");

    double T     = readValue<double>("  Truncamento T   (default 3.5) = ");
    double hStep = readValue<double>("  Passo      h   (default 0.05) = ");
    if (T     <= 0) T     = 3.5;
    if (hStep <= 0) hStep = 0.05;

    SingularityIntegrator integrator(fp, m, T, hStep);
    double result = integrator.compute(a, b);
    printResult(integrator.name(),
                "∫[" + std::to_string(a) + "," + std::to_string(b) + "] f(x)dx",
                result);
}

static void menuIntegracao() {
    while (true) {
        header("INTEGRAÇÃO NUMÉRICA");
        std::cout << "  [1] Newton-Cotes Fechado (Trapézio, Simpson)\n";
        std::cout << "  [2] Newton-Cotes Aberto  (Trapézio Aberto, Milne)\n";
        std::cout << "  [3] Gauss-Legendre       (n = 2..5)\n";
        std::cout << "  [4] Gauss Especial       (Hermite / Laguerre / Chebyshev)\n";
        std::cout << "  [5] Singularidades       (Exp. Simples / Dupla)\n";
        std::cout << "  [0] Voltar\n";

        int op = readInt("\n  Opção: ", 0, 5);
        if (op == 0) break;
        try {
            switch (op) {
                case 1: menuNewtonCotesFechado(); break;
                case 2: menuNewtonCotesAberto();  break;
                case 3: menuGaussLegendre();      break;
                case 4: menuGaussEspecial();      break;
                case 5: menuSingularidades();     break;
            }
        } catch (const std::exception& e) {
            std::cout << "\n  [ERRO] " << e.what() << "\n";
        }
        std::cout << "\n  Pressione ENTER para continuar...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// =============================================================================
// MÓDULO 3 — AUTOVALORES E AUTOVETORES
// =============================================================================

/// Lê dimensões e elementos de uma matriz simétrica.
static Matrix lerMatrizQuadrada(const std::string& label = "A") {
    int n = readInt("\n  Dimensão n (matriz " + label + " é n×n): ", 1, 20);
    std::cout << "\n  Matriz " << label << " (" << n << "×" << n << "):\n";
    return readMatrix(n, n);
}

static void imprimirEigenResult(const EigenResult& r, const std::string& method) {
    clearLine();
    std::cout << "  Método    : " << method << "\n";
    std::cout << "  Autovalor : "
              << std::fixed << std::setprecision(10) << r.eigenvalue << "\n";
    printVec("Autovetor ", r.eigenvector, 8);
    std::cout << "  Iterações : " << r.iterations << "\n";
    std::cout << "  Convergiu : " << (r.converged ? "Sim" : "NÃO") << "\n";
    clearLine();
}

static void menuPotenciaRegular() {
    subheader("Potência Regular — Autovalor Dominante");
    Matrix A = lerMatrizQuadrada();
    A.print("\n  Matriz A:", 12, 4);

    double tol    = readValue<double>("\n  Tolerância ε = ");
    int    maxIt  = readInt("  Máx. iterações = ", 1, 100000);

    PowerRegular solver;
    EigenResult  r = solver.solve(A, tol, maxIt);
    imprimirEigenResult(r, solver.name());
}

static void menuPotenciaInverso() {
    subheader("Potência Inverso — Menor Autovalor");
    Matrix A = lerMatrizQuadrada();
    A.print("\n  Matriz A:", 12, 4);

    double tol   = readValue<double>("\n  Tolerância ε = ");
    int    maxIt = readInt("  Máx. iterações = ", 1, 100000);

    PowerInverse solver;
    EigenResult  r = solver.solve(A, tol, maxIt);
    imprimirEigenResult(r, solver.name());
}

static void menuPotenciaDeslocamento() {
    subheader("Potência com Deslocamento — Autovalor próximo de μ");
    Matrix A = lerMatrizQuadrada();
    A.print("\n  Matriz A:", 12, 4);

    double mu    = readValue<double>("\n  Deslocamento μ = ");
    double tol   = readValue<double>("  Tolerância  ε = ");
    int    maxIt = readInt("  Máx. iterações = ", 1, 100000);

    PowerShifted solver(mu);
    EigenResult  r = solver.solve(A, tol, maxIt);
    imprimirEigenResult(r, solver.name());
}

static void menuHouseholder() {
    subheader("Transformação de Householder — Tridiagonalização");
    Matrix A = lerMatrizQuadrada();
    A.print("\n  Matriz A:", 12, 4);

    HouseholderMethod hh;
    auto res = hh.tridiagonalize(A);

    std::cout << "\n  Matriz Tridiagonal T:\n";
    res.T.print("", 12, 6);

    std::cout << "\n  Matriz de Transformação Q  (QᵀAQ = T):\n";
    res.Q.print("", 12, 6);

    // Verificação: imprime QᵀAQ para conferir
    Matrix check = res.Q.transpose() * A * res.Q;
    std::cout << "\n  Verificação QᵀAQ (deve ser ≈ T):\n";
    check.print("", 12, 6);
    clearLine();
}

static void menuJacobi() {
    subheader("Método de Jacobi — Todos os Autovalores/Autovetores");
    Matrix A = lerMatrizQuadrada();
    A.print("\n  Matriz A:", 12, 4);

    double tol   = readValue<double>("\n  Tolerância ε = ");
    int    maxIt = readInt("  Máx. iterações = ", 1, 200000);

    JacobiMethod solver;
    JacobiResult r = solver.solve(A, tol, maxIt);

    clearLine();
    std::cout << "  Método    : " << solver.name() << "\n";
    std::cout << "  Iterações : " << r.iterations << "\n";
    std::cout << "  Convergiu : " << (r.converged ? "Sim" : "NÃO") << "\n\n";

    int n = (int)r.eigenvalues.size();
    for (int i = 0; i < n; ++i) {
        std::cout << "  λ" << i + 1 << " = "
                  << std::fixed << std::setprecision(10) << r.eigenvalues[i] << "\n";
        std::cout << "  v" << i + 1 << " = [";
        for (int j = 0; j < n; ++j) {
            std::cout << std::fixed << std::setprecision(6) << r.eigenvectors(j, i);
            if (j + 1 < n) std::cout << ",  ";
        }
        std::cout << "]\n\n";
    }
    clearLine();
}

static void menuQR() {
    subheader("Método QR — Todos os Autovalores/Autovetores");
    std::cout << "  [1] QR com Householder (pré-tridiagonalização, recomendado)\n";
    std::cout << "  [2] QR direto (sem pré-tridiagonalização)\n";
    int useHH = readInt("\n  Escolha: ", 1, 2);

    Matrix A = lerMatrizQuadrada();
    A.print("\n  Matriz A:", 12, 4);

    double tol   = readValue<double>("\n  Tolerância ε = ");
    int    maxIt = readInt("  Máx. iterações = ", 1, 200000);

    QRMethod solver(useHH == 1);
    QRResult r = solver.solve(A, tol, maxIt);

    clearLine();
    std::cout << "  Método    : " << solver.name() << "\n";
    std::cout << "  Iterações : " << r.iterations << "\n";
    std::cout << "  Convergiu : " << (r.converged ? "Sim" : "NÃO") << "\n\n";

    int n = (int)r.eigenvalues.size();
    for (int i = 0; i < n; ++i) {
        std::cout << "  λ" << i + 1 << " = "
                  << std::fixed << std::setprecision(10) << r.eigenvalues[i] << "\n";
        std::cout << "  v" << i + 1 << " = [";
        for (int j = 0; j < n; ++j) {
            std::cout << std::fixed << std::setprecision(6) << r.eigenvectors(j, i);
            if (j + 1 < n) std::cout << ",  ";
        }
        std::cout << "]\n\n";
    }
    clearLine();
}

static void menuSVD() {
    subheader("SVD — Decomposição em Valores Singulares  A = UΣVᵀ");
    int m = readInt("\n  Número de linhas   m = ", 1, 20);
    int n = readInt("  Número de colunas  n = ", 1, 20);
    std::cout << "\n  Matriz A (" << m << "×" << n << "):\n";
    Matrix A = readMatrix(m, n);
    A.print("\n  Matriz A:", 12, 4);

    double tol = readValue<double>("\n  Tolerância ε para valor singular nulo = ");

    SVDSolver solver(tol);
    SVDResult r = solver.solve(A);

    clearLine();
    std::cout << "  Método : " << solver.name() << "\n";
    std::cout << "  Posto  : " << r.rank << "\n\n";

    std::cout << "  Valores Singulares (σᵢ = √λᵢ):\n";
    for (int i = 0; i < r.rank; ++i)
        std::cout << "    σ" << i + 1 << " = "
                  << std::fixed << std::setprecision(10) << r.singularValues[i] << "\n";

    std::cout << "\n  Matriz U (" << m << "×" << r.rank << ") — vetores esquerdos:\n";
    r.U.print("", 12, 6);

    std::cout << "\n  Matriz V (" << n << "×" << r.rank << ") — vetores direitos:\n";
    r.V.print("", 12, 6);

    // Verificação numérica: ||A - UΣVᵀ||_F
    Matrix Sigma(r.rank, r.rank, 0.0);
    for (int i = 0; i < r.rank; ++i) Sigma(i, i) = r.singularValues[i];
    Matrix reconstructed = r.U * Sigma * r.V.transpose();
    Matrix diff = A - reconstructed;
    std::cout << "\n  ||A - UΣVᵀ||_F (erro de reconstrução) = "
              << std::scientific << std::setprecision(4) << diff.frobeniusNorm() << "\n";
    clearLine();
}

static void menuAutovalores() {
    while (true) {
        header("AUTOVALORES E AUTOVETORES");
        std::cout << "  [1] Potência Regular      — autovalor dominante\n";
        std::cout << "  [2] Potência Inverso      — menor autovalor\n";
        std::cout << "  [3] Potência c/ Deslocam. — autovalor próximo de μ\n";
        std::cout << "  [4] Householder           — tridiagonalização\n";
        std::cout << "  [5] Jacobi                — todos autovalores/vetores\n";
        std::cout << "  [6] QR                    — todos autovalores/vetores\n";
        std::cout << "  [7] SVD                   — decomposição singular\n";
        std::cout << "  [0] Voltar\n";

        int op = readInt("\n  Opção: ", 0, 7);
        if (op == 0) break;
        try {
            switch (op) {
                case 1: menuPotenciaRegular();    break;
                case 2: menuPotenciaInverso();    break;
                case 3: menuPotenciaDeslocamento();break;
                case 4: menuHouseholder();        break;
                case 5: menuJacobi();             break;
                case 6: menuQR();                 break;
                case 7: menuSVD();                break;
            }
        } catch (const std::exception& e) {
            std::cout << "\n  [ERRO] " << e.what() << "\n";
        }
        std::cout << "\n  Pressione ENTER para continuar...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// =============================================================================
// MENU PRINCIPAL
// =============================================================================

static void splash() {
    std::cout << R"(
 ╔══════════════════════════════════════════════════════════════════╗
 ║          MÉTODOS NUMÉRICOS — Interface de Linha de Comando       ║
 ║  Diferenciação · Integração · Autovalores/Autovetores · SVD      ║
 ╚══════════════════════════════════════════════════════════════════╝

  Funções aceitam expressões como:
    sin(x), cos(x), exp(x), log(x), sqrt(x), x^2, x*y, etc.
)";
}

int main() {
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);

    splash();
    while (true) {
        header("MENU PRINCIPAL");
        std::cout << "  [1] Diferenciação Numérica\n";
        std::cout << "  [2] Integração Numérica\n";
        std::cout << "  [3] Autovalores e Autovetores\n";
        std::cout << "  [0] Sair\n";

        int op = readInt("\n  Opção: ", 0, 3);
        if (op == 0) {
            std::cout << "\n  Até mais!\n\n";
            break;
        }
        switch (op) {
            case 1: menuDiferenciacao(); break;
            case 2: menuIntegracao();    break;
            case 3: menuAutovalores();   break;
        }
    }
    return 0;
}