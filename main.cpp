// =============================================================================
// main.cpp
// Interface de linha de comando (CLI) para explorar todos os métodos numéricos.
//
// Estrutura de menus:
//   [1] Diferenciação Numérica
//   [2] Integração Numérica
//   [3] Autovalores e Autovetores
//   [4] Problemas de Valor Inicial (PVI)
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

#ifdef _WIN32
#include <windows.h>
#endif

// Core
#include "core/FunctionParser.hpp"
#include "core/SystemParser.hpp"

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

// Problemas de Valor Inicial (PVI)
#include "ivp/EulerMethod.hpp"
#include "ivp/RungeKutta.hpp"
#include "ivp/AdamsBashforth.hpp"

// Utils
#include "Matrix.hpp"

using namespace nm;

// =============================================================================
// Helpers de UI
// =============================================================================

static void clearLine()
{
    std::cout << std::string(70, '-') << "\n";
}

static void header(const std::string &title)
{
    std::cout << "\n"
              << std::string(70, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(70, '=') << "\n";
}

static void subheader(const std::string &title)
{
    std::cout << "\n"
              << std::string(50, '-') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(50, '-') << "\n";
}

template <typename T>
T readValue(const std::string &prompt)
{
    T value{};
    while (true)
    {
        std::cout << prompt;
        std::string line;
        std::getline(std::cin, line);
        std::istringstream iss(line);
        if (iss >> value)
            break;
        std::cout << "  [!] Entrada inválida. Tente novamente.\n";
    }
    return value;
}

static std::string readString(const std::string &prompt)
{
    std::cout << prompt;
    std::string line;
    std::getline(std::cin, line);
    return line;
}

static int readInt(const std::string &prompt, int lo, int hi)
{
    while (true)
    {
        int v = readValue<int>(prompt);
        if (v >= lo && v <= hi)
            return v;
        std::cout << "  [!] Valor fora do intervalo [" << lo << ", " << hi << "].\n";
    }
}

static std::shared_ptr<FunctionParser> readFunction1D()
{
    auto fp = std::make_shared<FunctionParser>();
    while (true)
    {
        std::cout << "  f(x) = ";
        std::string expr;
        std::getline(std::cin, expr);
        try
        {
            fp->parse(expr);
            std::cout << "  → Função aceita: f(x) = " << expr << "\n";
            return fp;
        }
        catch (const std::exception &e)
        {
            std::cout << "  [!] " << e.what() << "\n";
        }
    }
}

static Matrix readMatrix(int rows, int cols)
{
    Matrix M(rows, cols);
    std::cout << "  Digite os elementos linha a linha (separados por espaço):\n";
    for (int i = 0; i < rows; ++i)
    {
        while (true)
        {
            std::cout << "  Linha " << i + 1 << ": ";
            std::string line;
            std::getline(std::cin, line);
            std::istringstream iss(line);
            bool ok = true;
            for (int j = 0; j < cols; ++j)
            {
                double v;
                if (!(iss >> v))
                {
                    ok = false;
                    break;
                }
                M(i, j) = v;
            }
            if (ok)
                break;
            std::cout << "  [!] Esperados " << cols << " valores. Tente novamente.\n";
        }
    }
    return M;
}

static void printVec(const std::string &label, const std::vector<double> &v, int prec = 8)
{
    std::cout << "  " << label << " = [";
    for (size_t i = 0; i < v.size(); ++i)
    {
        std::cout << std::fixed << std::setprecision(prec) << v[i];
        if (i + 1 < v.size())
            std::cout << ",  ";
    }
    std::cout << "]\n";
}

static void printResult(const std::string &method, const std::string &label, double value, int prec = 10)
{
    clearLine();
    std::cout << "  Método  : " << method << "\n";
    std::cout << "  " << label << " = " << std::fixed << std::setprecision(prec) << value << "\n";
    clearLine();
}

// =============================================================================
// MÓDULO 1 — DIFERENCIAÇÃO NUMÉRICA
// =============================================================================

static void menuPrimeiraDerivada()
{
    subheader("Derivada Primeira f'(x)");
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
        FirstDerivMethod::FORWARD_2PT, FirstDerivMethod::FORWARD_3PT,
        FirstDerivMethod::FORWARD_4PT, FirstDerivMethod::BACKWARD_2PT,
        FirstDerivMethod::BACKWARD_3PT, FirstDerivMethod::BACKWARD_4PT,
        FirstDerivMethod::CENTRAL_2PT, FirstDerivMethod::CENTRAL_4PT};
    FirstDerivMethod m = methods[choice - 1];

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();

    double x = readValue<double>("\n  Ponto x = ");
    double h = readValue<double>("  Passo  h = ");

    FirstDerivative deriv(fp, m);
    double result = deriv.compute(x, h);
    printResult(deriv.name(), "f'(" + std::to_string(x) + ")", result);
}

static void menuSegundaDerivada()
{
    subheader("Derivada Segunda f''(x)");
    std::cout << "  [1] Central 3 pontos  O(h²)\n";
    std::cout << "  [2] Central 5 pontos  O(h⁴)\n";

    int choice = readInt("\n  Escolha o método: ", 1, 2);
    SecondDerivMethod m = (choice == 1) ? SecondDerivMethod::CENTRAL_3PT
                                        : SecondDerivMethod::CENTRAL_5PT;

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double x = readValue<double>("\n  Ponto x = ");
    double h = readValue<double>("  Passo  h = ");

    SecondDerivative deriv(fp, m);
    double result = deriv.compute(x, h);
    printResult(deriv.name(), "f''(" + std::to_string(x) + ")", result);
}

static void menuTerceiraDerivada()
{
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

static void menuDiferenciacao()
{
    while (true)
    {
        header("DIFERENCIAÇÃO NUMÉRICA");
        std::cout << "  [1] Derivada Primeira\n";
        std::cout << "  [2] Derivada Segunda\n";
        std::cout << "  [3] Derivada Terceira\n";
        std::cout << "  [0] Voltar\n";

        int op = readInt("\n  Opção: ", 0, 4);
        if (op == 0)
            break;
        try
        {
            switch (op)
            {
            case 1:
                menuPrimeiraDerivada();
                break;
            case 2:
                menuSegundaDerivada();
                break;
            case 3:
                menuTerceiraDerivada();
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cout << "\n  [ERRO] " << e.what() << "\n";
        }
        std::cout << "\n  Pressione ENTER para continuar...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// =============================================================================
// MÓDULO 2 — INTEGRAÇÃO NUMÉRICA
// =============================================================================

static void menuNewtonCotesFechado()
{
    subheader("Newton-Cotes Fechado");
    std::cout << "  [1] Regra do Trapézio\n";
    std::cout << "  [2] Simpson 1/3\n";
    std::cout << "  [3] Simpson 3/8\n";
    int choice = readInt("\n  Escolha: ", 1, 3);
    NewtonCotesMethod methods[] = {NewtonCotesMethod::TRAPEZOID, NewtonCotesMethod::SIMPSON_1_3, NewtonCotesMethod::SIMPSON_3_8};

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double a = readValue<double>("\n  Limite inferior a = ");
    double b = readValue<double>("  Limite superior b = ");

    NewtonCotes integrator(fp, methods[choice - 1]);
    double result = integrator.compute(a, b);
    printResult(integrator.name(), "∫[" + std::to_string(a) + "," + std::to_string(b) + "] f(x)dx", result);
}

static void menuNewtonCotesAberto()
{
    subheader("Newton-Cotes Aberto");
    std::cout << "  [1] Trapézio Aberto\n";
    std::cout << "  [2] Regra de Milne\n";
    int choice = readInt("\n  Escolha: ", 1, 2);
    NewtonCotesMethod m = (choice == 1) ? NewtonCotesMethod::OPEN_TRAPEZOID : NewtonCotesMethod::MILNE;

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double a = readValue<double>("\n  Limite inferior a = ");
    double b = readValue<double>("  Limite superior b = ");

    NewtonCotes integrator(fp, m);
    double result = integrator.compute(a, b);
    printResult(integrator.name(), "∫[" + std::to_string(a) + "," + std::to_string(b) + "] f(x)dx", result);
}

static void menuGaussLegendre()
{
    subheader("Quadratura de Gauss-Legendre");
    std::cout << "  Número de pontos n = 2, 3, 4 ou 5\n";
    int n = readInt("\n  n = ", 2, 5);

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double a = readValue<double>("\n  Limite inferior a = ");
    double b = readValue<double>("  Limite superior b = ");

    GaussLegendre integrator(fp, n);
    double result = integrator.compute(a, b);
    printResult(integrator.name(), "∫[" + std::to_string(a) + "," + std::to_string(b) + "] f(x)dx", result);
}

static void menuGaussEspecial()
{
    subheader("Quadraturas Especiais de Gauss (5 pontos)");
    std::cout << "  [1] Gauss-Hermite   — ∫_{-∞}^{+∞} e^{-x²} f(x) dx\n";
    std::cout << "  [2] Gauss-Laguerre  — ∫_{0}^{+∞}  e^{-x}  f(x) dx\n";
    std::cout << "  [3] Gauss-Chebyshev — ∫_{-1}^{1}  f(x)/√(1-x²) dx\n";
    int choice = readInt("\n  Escolha: ", 1, 3);
    SpecialGaussType types[] = {SpecialGaussType::HERMITE, SpecialGaussType::LAGUERRE, SpecialGaussType::CHEBYSHEV};

    int n = 5;
    if (choice == 3)
    {
        std::cout << "  Número de pontos de Chebyshev (1-20): ";
        n = readInt("", 1, 20);
    }

    std::cout << "\n  ATENÇÃO: digite apenas a parte f(x) (sem o peso).\n";
    std::cout << "\n  Digite f(x):\n";
    auto fp = readFunction1D();
    SpecialGauss integrator(fp, types[choice - 1], n);
    printResult(integrator.name(), "I", integrator.compute());
}

static void menuSingularidades()
{
    subheader("Estratégias para Singularidades");
    std::cout << "  [1] Exponencial Simples  [x(s) = m + r·tanh(s)]\n";
    std::cout << "  [2] Exponencial Dupla    [x(s) = m + r·tanh(π/2·sinh(s))]\n";
    int choice = readInt("\n  Escolha: ", 1, 2);
    SingularityMethod m = (choice == 1) ? SingularityMethod::SINGLE_EXP : SingularityMethod::DOUBLE_EXP;

    std::cout << "\n  Digite a função f(x):\n";
    auto fp = readFunction1D();
    double a = readValue<double>("\n  Limite inferior a = ");
    double b = readValue<double>("  Limite superior b = ");
    double T = readValue<double>("  Truncamento T   (default 3.5) = ");
    double hStep = readValue<double>("  Passo      h   (default 0.05) = ");
    if (T <= 0)
        T = 3.5;
    if (hStep <= 0)
        hStep = 0.05;

    SingularityIntegrator integrator(fp, m, T, hStep);
    printResult(integrator.name(), "∫[" + std::to_string(a) + "," + std::to_string(b) + "] f(x)dx", integrator.compute(a, b));
}

static void menuIntegracao()
{
    while (true)
    {
        header("INTEGRAÇÃO NUMÉRICA");
        std::cout << "  [1] Newton-Cotes Fechado (Trapézio, Simpson)\n";
        std::cout << "  [2] Newton-Cotes Aberto  (Trapézio Aberto, Milne)\n";
        std::cout << "  [3] Gauss-Legendre       (n = 2..5)\n";
        std::cout << "  [4] Gauss Especial       (Hermite / Laguerre / Chebyshev)\n";
        std::cout << "  [5] Singularidades       (Exp. Simples / Dupla)\n";
        std::cout << "  [0] Voltar\n";

        int op = readInt("\n  Opção: ", 0, 5);
        if (op == 0)
            break;
        try
        {
            switch (op)
            {
            case 1:
                menuNewtonCotesFechado();
                break;
            case 2:
                menuNewtonCotesAberto();
                break;
            case 3:
                menuGaussLegendre();
                break;
            case 4:
                menuGaussEspecial();
                break;
            case 5:
                menuSingularidades();
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cout << "\n  [ERRO] " << e.what() << "\n";
        }
        std::cout << "\n  Pressione ENTER para continuar...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// =============================================================================
// MÓDULO 3 — AUTOVALORES E AUTOVETORES
// =============================================================================

static Matrix lerMatrizQuadrada(const std::string &label = "A")
{
    int n = readInt("\n  Dimensão n (matriz " + label + " é n×n): ", 1, 20);
    std::cout << "\n  Matriz " << label << " (" << n << "×" << n << "):\n";
    return readMatrix(n, n);
}

static void imprimirEigenResult(const EigenResult &r, const std::string &method)
{
    clearLine();
    std::cout << "  Método    : " << method << "\n";
    std::cout << "  Autovalor : " << std::fixed << std::setprecision(10) << r.eigenvalue << "\n";
    printVec("Autovetor ", r.eigenvector, 8);
    std::cout << "  Iterações : " << r.iterations << "\n";
    std::cout << "  Convergiu : " << (r.converged ? "Sim" : "NÃO") << "\n";
    clearLine();
}

static void menuPotenciaRegular()
{
    subheader("Potência Regular");
    Matrix A = lerMatrizQuadrada();
    double tol = readValue<double>("\n  Tolerância ε = ");
    int maxIt = readInt("  Máx. iterações = ", 1, 100000);
    PowerRegular solver;
    imprimirEigenResult(solver.solve(A, tol, maxIt), solver.name());
}

static void menuPotenciaInverso()
{
    subheader("Potência Inverso");
    Matrix A = lerMatrizQuadrada();
    double tol = readValue<double>("\n  Tolerância ε = ");
    int maxIt = readInt("  Máx. iterações = ", 1, 100000);
    PowerInverse solver;
    imprimirEigenResult(solver.solve(A, tol, maxIt), solver.name());
}

static void menuPotenciaDeslocamento()
{
    subheader("Potência com Deslocamento");
    Matrix A = lerMatrizQuadrada();
    double mu = readValue<double>("\n  Deslocamento μ = ");
    double tol = readValue<double>("  Tolerância  ε = ");
    int maxIt = readInt("  Máx. iterações = ", 1, 100000);
    PowerShifted solver(mu);
    imprimirEigenResult(solver.solve(A, tol, maxIt), solver.name());
}

static void menuHouseholder()
{
    subheader("Householder — Tridiagonalização");
    Matrix A = lerMatrizQuadrada();
    HouseholderMethod hh;
    auto res = hh.tridiagonalize(A);
    std::cout << "\n  Matriz Tridiagonal T:\n";
    res.T.print("", 12, 6);
    std::cout << "\n  Matriz Q:\n";
    res.Q.print("", 12, 6);
    clearLine();
}

static void menuJacobi()
{
    subheader("Método de Jacobi");
    Matrix A = lerMatrizQuadrada();
    double tol = readValue<double>("\n  Tolerância ε = ");
    int maxIt = readInt("  Máx. iterações = ", 1, 200000);
    JacobiMethod solver;
    JacobiResult r = solver.solve(A, tol, maxIt);

    clearLine();
    std::cout << "  Método    : " << solver.name() << "\n";
    std::cout << "  Iterações : " << r.iterations << "\n\n";
    int n = (int)r.eigenvalues.size();
    for (int i = 0; i < n; ++i)
    {
        std::cout << "  λ" << i + 1 << " = " << std::fixed << std::setprecision(10) << r.eigenvalues[i] << "\n";
        std::cout << "  v" << i + 1 << " = [";
        for (int j = 0; j < n; ++j)
        {
            std::cout << std::fixed << std::setprecision(6) << r.eigenvectors(j, i);
            if (j + 1 < n)
                std::cout << ",  ";
        }
        std::cout << "]\n\n";
    }
    clearLine();
}

static void menuQR()
{
    subheader("Método QR");
    std::cout << "  [1] QR com Householder\n  [2] QR direto\n";
    int useHH = readInt("\n  Escolha: ", 1, 2);
    Matrix A = lerMatrizQuadrada();
    double tol = readValue<double>("\n  Tolerância ε = ");
    int maxIt = readInt("  Máx. iterações = ", 1, 200000);

    QRMethod solver(useHH == 1);
    QRResult r = solver.solve(A, tol, maxIt);

    clearLine();
    std::cout << "  Método    : " << solver.name() << "\n";
    std::cout << "  Iterações : " << r.iterations << "\n\n";
    int n = (int)r.eigenvalues.size();
    for (int i = 0; i < n; ++i)
    {
        std::cout << "  λ" << i + 1 << " = " << std::fixed << std::setprecision(10) << r.eigenvalues[i] << "\n";
        std::cout << "  v" << i + 1 << " = [";
        for (int j = 0; j < n; ++j)
        {
            std::cout << std::fixed << std::setprecision(6) << r.eigenvectors(j, i);
            if (j + 1 < n)
                std::cout << ",  ";
        }
        std::cout << "]\n\n";
    }
    clearLine();
}

static void menuSVD()
{
    subheader("SVD");
    int m = readInt("\n  Linhas m = ", 1, 20);
    int n = readInt("  Colunas n = ", 1, 20);
    Matrix A = readMatrix(m, n);
    double tol = readValue<double>("\n  Tolerância ε = ");

    SVDSolver solver(tol);
    SVDResult r = solver.solve(A);

    clearLine();
    std::cout << "  Método : " << solver.name() << "\n  Posto  : " << r.rank << "\n\n";
    for (int i = 0; i < r.rank; ++i)
        std::cout << "    σ" << i + 1 << " = " << std::fixed << std::setprecision(10) << r.singularValues[i] << "\n";
    std::cout << "\n  Matriz U:\n";
    r.U.print("", 12, 6);
    std::cout << "\n  Matriz V:\n";
    r.V.print("", 12, 6);
    clearLine();
}

static void menuAutovalores()
{
    while (true)
    {
        header("AUTOVALORES E AUTOVETORES");
        std::cout << "  [1] Potência Regular\n  [2] Potência Inverso\n  [3] Potência c/ Deslocamento\n";
        std::cout << "  [4] Householder\n  [5] Jacobi\n  [6] QR\n  [7] SVD\n  [0] Voltar\n";
        int op = readInt("\n  Opção: ", 0, 7);
        if (op == 0)
            break;
        try
        {
            switch (op)
            {
            case 1:
                menuPotenciaRegular();
                break;
            case 2:
                menuPotenciaInverso();
                break;
            case 3:
                menuPotenciaDeslocamento();
                break;
            case 4:
                menuHouseholder();
                break;
            case 5:
                menuJacobi();
                break;
            case 6:
                menuQR();
                break;
            case 7:
                menuSVD();
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cout << "\n  [ERRO] " << e.what() << "\n";
        }
        std::cout << "\n  Pressione ENTER para continuar...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// =============================================================================
// MÓDULO 4 — PROBLEMAS DE VALOR INICIAL (PVI)
// =============================================================================

// Helper local para formatar e imprimir os resultados do PVI de forma limpa.
static void printIVPResult(const nm::IVPResult &res, const std::string &solverName)
{
    clearLine();
    std::cout << "  Solver    : " << solverName << "\n";
    std::cout << "  Passos    : " << res.steps << "\n";
    std::cout << "  Status    : " << (res.completed ? "Concluído com sucesso" : "Interrompido/Limite atingido") << "\n";

    int n = res.y.empty() ? 0 : (int)res.y[0].size();

    // Imprimir sempre o estado final
    std::cout << "\n  [Estado Final em t = " << std::fixed << std::setprecision(6) << res.t.back() << "]\n";
    for (int i = 0; i < n; ++i)
    {
        std::cout << "    s" << i << " = " << std::setprecision(10) << res.y.back()[i] << "\n";
    }
    clearLine();

    // Se a tabela for pequena, exibimos na íntegra. Caso contrário, evitamos poluir o terminal.
    if (res.steps <= 50 && res.steps > 0)
    {
        std::cout << "\n  Tabela Completa (t, s0, s1, ...):\n";
        for (size_t k = 0; k < res.t.size(); ++k)
        {
            std::cout << "  t=" << std::setw(8) << std::left << std::setprecision(4) << res.t[k] << " | ";
            for (int i = 0; i < n; ++i)
            {
                std::cout << std::setw(12) << std::right << std::setprecision(6) << res.y[k][i] << " ";
            }
            std::cout << "\n";
        }
        clearLine();
    }
    else if (res.steps > 50)
    {
        std::cout << "\n  (Tabela completa com " << res.steps << " linhas foi omitida para limpeza do terminal.)\n";
        std::cout << "  (Apenas o estado final acima foi exibido.)\n";
        clearLine();
    }
}

static void menuPVI()
{
    while (true)
    {
        header("PROBLEMAS DE VALOR INICIAL (PVI)");
        std::cout << "  [1] Euler Explícito\n";
        std::cout << "  [2] Euler Implícito (Ponto Fixo)\n";
        std::cout << "  [3] Runge-Kutta 2ª Ordem (Heun)\n";
        std::cout << "  [4] Runge-Kutta 3ª Ordem\n";
        std::cout << "  [5] Runge-Kutta 4ª Ordem (Clássico)\n";
        std::cout << "  [6] Adams-Bashforth-Moulton 2ª Ordem\n";
        std::cout << "  [7] Adams-Bashforth-Moulton 3ª Ordem\n";
        std::cout << "  [8] Adams-Bashforth-Moulton 4ª Ordem\n";
        std::cout << "  [0] Voltar\n";

        int op = readInt("\n  Opção (Escolha o Solucionador): ", 0, 8);
        if (op == 0)
            break;

        try
        {
            subheader("Configuração do Sistema de EDOs");
            int n = readInt("\n  Dimensão do sistema (n): ", 1, 20);

            std::vector<std::string> exprs(n);
            std::cout << "\n  Digite as expressões usando as variáveis de estado (s0, s1...) e tempo (t).\n";
            for (int i = 0; i < n; ++i)
            {
                exprs[i] = readString("  F_" + std::to_string(i) + "(s, t) = ");
            }

            // Inicializa o SystemParser
            nm::SystemParser parser;
            parser.setup(n, exprs);
            std::cout << "  → Sistema processado com sucesso.\n";

            subheader("Condições Iniciais e Parâmetros");
            std::vector<double> s0(n);
            for (int i = 0; i < n; ++i)
            {
                s0[i] = readValue<double>("  Condição inicial s" + std::to_string(i) + "(t0) = ");
            }

            double t0 = readValue<double>("  Tempo inicial (t0) = ");
            double tf = readValue<double>("  Tempo final (tf) = ");
            double h = readValue<double>("  Tamanho do passo (h) = ");

            // Seleção de Solver usando polimorfismo
            std::unique_ptr<nm::IVPSolver> solver;
            switch (op)
            {
            case 1:
                solver = std::make_unique<nm::EulerExplicit>();
                break;
            case 2:
                solver = std::make_unique<nm::EulerImplicit>();
                break; // Aceita tolerâncias padrão
            case 3:
                solver = std::make_unique<nm::RungeKutta2>();
                break;
            case 4:
                solver = std::make_unique<nm::RungeKutta3>();
                break;
            case 5:
                solver = std::make_unique<nm::RungeKutta4>();
                break;
            case 6:
                solver = std::make_unique<nm::AdamsBashforth2>();
                break;
            case 7:
                solver = std::make_unique<nm::AdamsBashforth3>();
                break;
            case 8:
                solver = std::make_unique<nm::AdamsBashforth4>();
                break;
            }

            std::cout << "\n  [Executando Solver...]\n";
            nm::IVPResult result = solver->solve(parser, s0, t0, tf, h);

            // Impressão robusta
            printIVPResult(result, solver->name());
        }
        catch (const std::exception &e)
        {
            std::cout << "\n  [ERRO] " << e.what() << "\n";
        }

        std::cout << "\n  Pressione ENTER para continuar...";
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

// =============================================================================
// MENU PRINCIPAL
// =============================================================================

static void splash()
{
    std::cout << R"(
 ╔══════════════════════════════════════════════════════════════════╗
 ║          MÉTODOS NUMÉRICOS — Interface de Linha de Comando       ║
 ║ Diferenciação · Integração · Autovalores/Vetores · SVD · PVI (EDO)║
 ╚══════════════════════════════════════════════════════════════════╝
)";
}

int main()
{
#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif

    splash();
    while (true)
    {
        header("MENU PRINCIPAL");
        std::cout << "  [1] Diferenciação Numérica\n";
        std::cout << "  [2] Integração Numérica\n";
        std::cout << "  [3] Autovalores e Autovetores\n";
        std::cout << "  [4] Problemas de Valor Inicial (PVI)\n";
        std::cout << "  [0] Sair\n";

        int op = readInt("\n  Opção: ", 0, 4);
        if (op == 0)
        {
            std::cout << "\n  Até mais!\n\n";
            break;
        }
        switch (op)
        {
        case 1:
            menuDiferenciacao();
            break;
        case 2:
            menuIntegracao();
            break;
        case 3:
            menuAutovalores();
            break;
        case 4:
            menuPVI();
            break;
        }
    }
    return 0;
}