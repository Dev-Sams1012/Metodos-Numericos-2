#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <functional>

#ifdef _WIN32
#include <windows.h>
#endif

#include "core/FunctionParser.hpp"
#include "core/SystemParser.hpp"

#include "differentiation/FirstDerivative.hpp"
#include "differentiation/SecondDerivative.hpp"
#include "differentiation/ThirdDerivative.hpp"

#include "integration/NewtonCotes.hpp"
#include "integration/GaussLegendre.hpp"
#include "integration/SpecialGauss.hpp"
#include "integration/SingularityStrategies.hpp"

#include "eigenvalues/PowerMethod.hpp"
#include "eigenvalues/HouseholderMethod.hpp"
#include "eigenvalues/JacobiMethod.hpp"
#include "eigenvalues/QRMethod.hpp"
#include "eigenvalues/SVDMethod.hpp"

#include "ivp/EulerMethod.hpp"
#include "ivp/RungeKutta.hpp"
#include "ivp/AdamsBashforth.hpp"

#include "bvp/FDMethod1D.hpp"
#include "bvp/FDMethod2D.hpp"

#include "Matrix.hpp"

using namespace nm;

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
    while (true)
    {
        std::cout << prompt;

        std::string expr;
        std::getline(std::cin, expr);

        try
        {
            FunctionParser fp;
            fp.parse(expr);
            return fp.eval(0.0);
        }
        catch (...)
        {
            std::cout << "  [!] Entrada inválida. Tente novamente.\n";
        }
    }
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

static void menuNewtonCotesFechado()
{
    subheader("Newton-Cotes Fechado");
    std::cout << "  [1] Regra do Trapézio (grau 1)\n";
    std::cout << "  [2] Simpson 1/3 (grau 2)\n";
    std::cout << "  [3] Simpson 3/8 (grau 3)\n";
    std::cout << "  [4] Regra de Boole (grau 4)\n";
    int choice = readInt("\n  Escolha: ", 1, 4);
    NewtonCotesMethod methods[] = {NewtonCotesMethod::TRAPEZOID, NewtonCotesMethod::SIMPSON_1_3, NewtonCotesMethod::SIMPSON_3_8, NewtonCotesMethod::BOOLE};

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
    std::cout << "  [1] Trapézio Aberto (grau 1)\n";
    std::cout << "  [2] Regra de Milne (grau 2)\n";
    std::cout << "  [3] Aberta grau 3 (4 pontos)\n";
    std::cout << "  [4] Aberta grau 4 (5 pontos)\n";
    int choice = readInt("\n  Escolha: ", 1, 4);
    NewtonCotesMethod methodsOpen[] = {NewtonCotesMethod::OPEN_TRAPEZOID, NewtonCotesMethod::MILNE, NewtonCotesMethod::OPEN_DEGREE3, NewtonCotesMethod::OPEN_DEGREE4};
    NewtonCotesMethod m = methodsOpen[choice - 1];

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
    subheader("Quadraturas Especiais de Gauss");
    std::cout << "  [1] Gauss-Hermite   — ∫_{-∞}^{+∞} e^{-x²} f(x) dx\n";
    std::cout << "  [2] Gauss-Laguerre  — ∫_{0}^{+∞}  e^{-x}  f(x) dx\n";
    std::cout << "  [3] Gauss-Chebyshev — ∫_{-1}^{1}  f(x)/√(1-x²) dx\n";
    int choice = readInt("\n  Escolha: ", 1, 3);
    SpecialGaussType types[] = {SpecialGaussType::HERMITE, SpecialGaussType::LAGUERRE, SpecialGaussType::CHEBYSHEV};

    int n = 5;
    if (choice == 1 || choice == 2)
    {
        std::cout << "  Número de pontos n = 2, 3, 4 ou 5\n";
        n = readInt("  n = ", 2, 5);
    }
    else
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
        std::cout << "  [1] Newton-Cotes Fechado (graus 1 a 4)\n";
        std::cout << "  [2] Newton-Cotes Aberto  (graus 1 a 4)\n";
        std::cout << "  [3] Gauss-Legendre       (n = 2..5)\n";
        std::cout << "  [4] Gauss Especial       (Hermite / Laguerre n=2..5, Chebyshev n=1..20)\n";
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

    if (r.rank > 0)
    {
        Matrix Sigma(r.rank, r.rank, 0.0);
        for (int i = 0; i < r.rank; ++i)
            Sigma(i, i) = r.singularValues[i];

        Matrix Reconstructed = r.U * Sigma * r.V.transpose();

        std::cout << "\n  Matriz Original M:\n";
        A.print("", 12, 6);
        std::cout << "\n  Matriz Reconstruída (U * Σ * Vᵀ):\n";
        Reconstructed.print("", 12, 6);

        double maxError = 0.0;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
            {
                double err = std::abs(A(i, j) - Reconstructed(i, j));
                if (err > maxError)
                    maxError = err;
            }

        std::cout << "\n  Erro Máximo Absoluto (|M - UΣVᵀ|): "
                  << std::scientific << std::setprecision(4) << maxError << "\n";
    }

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

static void printIVPResult(const nm::IVPResult &res, const std::string &solverName)
{
    clearLine();
    std::cout << "  Solver    : " << solverName << "\n";
    std::cout << "  Passos    : " << res.steps << "\n";
    std::cout << "  Status    : " << (res.completed ? "Concluído com sucesso" : "Interrompido/Limite atingido") << "\n";

    int n = res.y.empty() ? 0 : (int)res.y[0].size();

    std::cout << "\n  [Estado Final em t = " << std::fixed << std::setprecision(6) << res.t.back() << "]\n";
    for (int i = 0; i < n; ++i)
        std::cout << "    s" << i << " = " << std::setprecision(10) << res.y.back()[i] << "\n";
    clearLine();

    if (res.steps <= 50 && res.steps > 0)
    {
        std::cout << "\n  Tabela Completa (t, s0, s1, ...):\n";
        for (size_t k = 0; k < res.t.size(); ++k)
        {
            std::cout << "  t=" << std::setw(8) << std::left << std::setprecision(4) << res.t[k] << " | ";
            for (int i = 0; i < n; ++i)
                std::cout << std::setw(12) << std::right << std::setprecision(6) << res.y[k][i] << " ";
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
                exprs[i] = readString("  F_" + std::to_string(i) + "(s, t) = ");

            nm::SystemParser parser;
            parser.setup(n, exprs);
            std::cout << "  → Sistema processado com sucesso.\n";

            subheader("Condições Iniciais e Parâmetros");
            std::vector<double> s0(n);
            for (int i = 0; i < n; ++i)
                s0[i] = readValue<double>("  Condição inicial s" + std::to_string(i) + "(t0) = ");

            double t0 = readValue<double>("  Tempo inicial (t0) = ");
            double tf = readValue<double>("  Tempo final (tf) = ");
            double h = readValue<double>("  Tamanho do passo (h) = ");

            std::unique_ptr<nm::IVPSolver> solver;
            switch (op)
            {
            case 1:
                solver = std::make_unique<nm::EulerExplicit>();
                break;
            case 2:
                solver = std::make_unique<nm::EulerImplicit>();
                break;
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

static void printBVPResult1D(const nm::BVPResult &r, const std::string &solverName)
{
    clearLine();
    std::cout << "  Solver   : " << solverName << "\n";
    std::cout << "  Divisões : n = " << r.n
              << "   (h = " << std::fixed << std::setprecision(6)
              << (r.x.back() - r.x.front()) / r.n << ")\n";
    std::cout << "  Status   : " << (r.solved ? "Sistema resolvido com sucesso" : "FALHOU") << "\n\n";

    if (!r.solved)
    {
        clearLine();
        return;
    }

    std::cout << "  " << std::setw(5) << "i"
              << std::setw(12) << "xi"
              << std::setw(20) << "u_aprox(xi)" << "\n";
    std::cout << "  " << std::string(37, '-') << "\n";

    for (int i = 0; i <= r.n; ++i)
        std::cout << "  " << std::setw(5) << i
                  << std::setw(12) << std::fixed << std::setprecision(6) << r.x[i]
                  << std::setw(20) << std::fixed << std::setprecision(10) << r.u[i]
                  << "\n";
    clearLine();
}

static void printBVPResult2D(const nm::BVPResult &r,
                             const std::string &solverName,
                             int nx, int ny)
{
    clearLine();
    std::cout << "  Solver   : " << solverName << "\n";
    std::cout << "  Grade    : " << nx << " × " << ny
              << "  →  " << (nx - 1) * (ny - 1) << " nós interiores\n";
    std::cout << "  Status   : " << (r.solved ? "Sistema resolvido com sucesso" : "FALHOU") << "\n\n";

    if (!r.solved)
    {
        clearLine();
        return;
    }

    const int mx = nx - 1;
    const int my = ny - 1;

    std::cout << "  u(x,y)      ";
    for (int i = 1; i <= mx; ++i)
        std::cout << std::setw(11) << std::fixed << std::setprecision(4) << r.x[i - 1];
    std::cout << "   ← x\n";
    std::cout << "  " << std::string(12 + 11 * mx, '-') << "\n";

    for (int j = my; j >= 1; --j)
    {
        double yj = r.y[(j - 1) * mx];
        std::cout << "  y=" << std::fixed << std::setprecision(4) << yj << " |";
        for (int i = 1; i <= mx; ++i)
        {
            int k = (j - 1) * mx + (i - 1);
            std::cout << std::setw(11) << std::fixed << std::setprecision(6) << r.u[k];
        }
        std::cout << "\n";
    }
    clearLine();
}

static nm::BoundaryCondition readBC(const std::string &side)
{
    std::cout << "\n  Condição em " << side << ":\n";
    std::cout << "    [1] Dirichlet  u(" << side << ") = valor\n";
    std::cout << "    [2] Neumann    u'(" << side << ") = valor\n";
    int tipo = readInt("  Tipo: ", 1, 2);
    double val = readValue<double>("  Valor: ");
    return (tipo == 1) ? nm::BoundaryCondition::dirichlet(val)
                       : nm::BoundaryCondition::neumann(val);
}

static void menuBVP1D()
{
    subheader("PVC 1D — Diferenças Finitas (CCs mistas D/N)");
    std::cout << "\n  Resolve: α·u''(x) + β·u'(x) + γ·u(x) = f(x)\n";
    std::cout << "  Condição em cada extremidade: Dirichlet u=valor  ou  Neumann u'=valor\n\n";

    double alpha = readValue<double>("  Coeficiente α (de u'') = ");
    double beta = readValue<double>("  Coeficiente β (de u')  = ");
    double gamma = readValue<double>("  Coeficiente γ (de u)   = ");

    std::cout << "\n  Função f(x) [ENTER para f ≡ 0]: ";
    std::string expr;
    std::getline(std::cin, expr);
    std::shared_ptr<nm::FunctionParser> fp;
    if (!expr.empty())
    {
        fp = std::make_shared<nm::FunctionParser>();
        try
        {
            fp->parse(expr);
            std::cout << "  → f(x) = " << expr << "\n";
        }
        catch (const std::exception &e)
        {
            std::cout << "  [!] " << e.what() << " — usando f ≡ 0.\n";
            fp.reset();
        }
    }

    double a = readValue<double>("\n  a (extremo esquerdo) = ");
    double b = readValue<double>("  b (extremo direito)  = ");

    nm::BoundaryCondition bcA = readBC("x=a");
    nm::BoundaryCondition bcB = readBC("x=b");

    int n = readInt("\n  Divisões n (≥ 2): ", 2, 10000);

    nm::FDMethod1D solver(alpha, beta, gamma);
    nm::BVPResult res = solver.solve(fp, a, b, bcA, bcB, n);

    auto bcLabel = [](const nm::BoundaryCondition &bc, const std::string &x)
    {
        return (bc.type == nm::BCType::DIRICHLET)
                   ? "u(" + x + ")=" + std::to_string(bc.value)
                   : "u'(" + x + ")=" + std::to_string(bc.value);
    };
    clearLine();
    std::cout << "  Solver   : " << solver.name() << "\n";
    std::cout << "  CC esq.  : " << bcLabel(bcA, "a") << "\n";
    std::cout << "  CC dir.  : " << bcLabel(bcB, "b") << "\n";
    std::cout << "  Divisões : n = " << res.n
              << "   (h = " << std::fixed << std::setprecision(6)
              << (res.x.back() - res.x.front()) / res.n << ")\n";
    std::cout << "  Status   : "
              << (res.solved ? "Sistema resolvido com sucesso" : "FALHOU") << "\n\n";

    if (res.solved)
    {
        std::cout << "  " << std::setw(5) << "i"
                  << std::setw(12) << "xi"
                  << std::setw(20) << "u_aprox(xi)" << "\n";
        std::cout << "  " << std::string(37, '-') << "\n";
        for (int i = 0; i <= res.n; ++i)
            std::cout << "  " << std::setw(5) << i
                      << std::setw(12) << std::fixed << std::setprecision(6) << res.x[i]
                      << std::setw(20) << std::fixed << std::setprecision(10) << res.u[i]
                      << "\n";
    }
    clearLine();

    std::cout << "  Comparar com solução exata? [s/n]: ";
    std::string resp;
    std::getline(std::cin, resp);
    if (resp.size() && (resp[0] == 's' || resp[0] == 'S'))
    {
        std::cout << "  u_exata(x) = ";
        std::string exprE;
        std::getline(std::cin, exprE);
        auto fpE = std::make_shared<nm::FunctionParser>();
        try
        {
            fpE->parse(exprE);
            clearLine();
            std::cout << "  " << std::setw(5) << "i"
                      << std::setw(12) << "xi"
                      << std::setw(18) << "u_aprox"
                      << std::setw(18) << "u_exata"
                      << std::setw(14) << "err_rel(%)" << "\n";
            std::cout << "  " << std::string(67, '-') << "\n";
            for (int i = 0; i <= res.n; ++i)
            {
                double xi = res.x[i];
                double ua_i = res.u[i];
                double ue_i = fpE->eval(xi);
                double err = (std::abs(ue_i) > 1e-14)
                                 ? std::abs(ua_i - ue_i) / std::abs(ue_i) * 100.0
                                 : std::abs(ua_i) * 100.0;
                std::cout << "  " << std::setw(5) << i
                          << std::setw(12) << std::fixed << std::setprecision(6) << xi
                          << std::setw(18) << std::setprecision(8) << ua_i
                          << std::setw(18) << std::setprecision(8) << ue_i
                          << std::setw(14) << std::setprecision(4) << err << "\n";
            }
            clearLine();
        }
        catch (const std::exception &e)
        {
            std::cout << "  [!] " << e.what() << "\n";
        }
    }
}

static void menuBVP2D()
{
    subheader("PVC 2D — Diferenças Finitas (Equação de Poisson)");
    std::cout << "\n  Resolve: ∂²u/∂x² + ∂²u/∂y² = f(x,y)\n";
    std::cout << "           u = 0 nas quatro bordas de [xa,xb]×[ya,yb]\n\n";
    std::cout << "  Exemplo da Aula 27 — PVC2:\n";
    std::cout << "    f(x,y) = 4,  domínio [0,1]×[0,1],  n=4\n\n";

    double xa = readValue<double>("  xa = ");
    double xb = readValue<double>("  xb = ");
    double ya = readValue<double>("  ya = ");
    double yb = readValue<double>("  yb = ");

    std::cout << "\n  f(x,y) — informe valor constante ou expressão em x\n";
    std::cout << "  [ENTER para f ≡ 0]: ";
    std::string fExpr;
    std::getline(std::cin, fExpr);

    std::function<double(double, double)> funcF;
    if (fExpr.empty())
    {
        funcF = [](double, double)
        { return 0.0; };
    }
    else
    {

        bool isConst = false;
        double fConst = 0.0;
        try
        {
            fConst = std::stod(fExpr);
            isConst = true;
        }
        catch (...)
        {
        }

        if (isConst)
        {
            funcF = [fConst](double, double)
            { return fConst; };
            std::cout << "  → f(x,y) = " << fConst << " (constante)\n";
        }
        else
        {
            auto fpF = std::make_shared<nm::FunctionParser>();
            try
            {
                fpF->parse(fExpr);
                funcF = [fpF](double x, double)
                { return fpF->eval(x); };
                std::cout << "  → f(x,y) = " << fExpr << " (expressão em x)\n";
            }
            catch (const std::exception &e)
            {
                std::cout << "  [!] " << e.what() << " — usando f ≡ 0.\n";
                funcF = [](double, double)
                { return 0.0; };
            }
        }
    }

    int nx = readInt("\n  Divisões em x (nx ≥ 2): ", 2, 200);
    int ny = readInt("  Divisões em y (ny, 0 = igual a nx): ", 0, 200);
    if (ny == 0)
        ny = nx;

    nm::FDMethod2D solver;
    nm::BVPResult res = solver.solve(funcF, xa, xb, ya, yb, nx, ny);
    printBVPResult2D(res, solver.name(), nx, ny);
}

static void menuBVP()
{
    while (true)
    {
        header("PROBLEMAS DE VALOR DE CONTORNO (PVC) — DIFERENÇAS FINITAS");
        std::cout << "  Método: grade uniforme + diferenças centrais O(h²) + sistema LU\n\n";
        std::cout << "  [1] PVC 1D  —  α·u'' + β·u' + γ·u = f(x),  u(a)=ua, u(b)=ub\n";
        std::cout << "  [2] PVC 2D  —  ∂²u/∂x² + ∂²u/∂y² = f(x,y)  (Poisson/Laplace)\n";
        std::cout << "  [0] Voltar\n";

        int op = readInt("\n  Opção: ", 0, 2);
        if (op == 0)
            break;

        try
        {
            switch (op)
            {
            case 1:
                menuBVP1D();
                break;
            case 2:
                menuBVP2D();
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

static void splash()
{
    std::cout << R"(
 ╔══════════════════════════════════════════════════════════════════╗
 ║          MÉTODOS NUMÉRICOS — Interface de Linha de Comando       ║
 ║  Diferenciação · Integração · Autovalores · PVI (EDO) · PVC(DF) ║
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
        std::cout << "  [4] Problemas de Valor Inicial  (PVI)\n";
        std::cout << "  [5] Problemas de Valor de Contorno (PVC) — Diferenças Finitas\n";
        std::cout << "  [0] Sair\n";

        int op = readInt("\n  Opção: ", 0, 5);
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
        case 5:
            menuBVP();
            break;
        }
    }
    return 0;
}