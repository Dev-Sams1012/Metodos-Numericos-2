

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

struct Image
{
    int width = 0, height = 0;
    std::vector<double> pixels;

    double &at(int i, int j)
    {
        return pixels[i * width + j];
    }

    double at(int i, int j) const
    {
        return pixels[i * width + j];
    }
};

Image readRealImage(const std::string &path)
{
    int width, height, channels;
    unsigned char *data = stbi_load(path.c_str(), &width, &height, &channels, 3);
    if (!data)
    {
        throw std::runtime_error("Não foi possível abrir/carregar a imagem: " + path);
    }

    Image img;
    img.width = width;
    img.height = height;
    img.pixels.resize(static_cast<size_t>(width) * height);

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            int idx = (i * width + j) * 3;
            double r = data[idx];
            double g = data[idx + 1];
            double b = data[idx + 2];

            img.at(i, j) = 0.299 * r + 0.587 * g + 0.114 * b;
        }
    }

    stbi_image_free(data);
    return img;
}

void writeRealImage(const std::string &path, const Image &img, bool binary = false)
{
    std::vector<unsigned char> out_pixels(static_cast<size_t>(img.width) * img.height);

    for (int i = 0; i < img.height; ++i)
    {
        for (int j = 0; j < img.width; ++j)
        {
            double v = img.at(i, j);
            if (binary)
            {
                v *= 255.0;
            }

            out_pixels[i * img.width + j] =
                static_cast<unsigned char>(std::clamp(v, 0.0, 255.0));
        }
    }

    int success = stbi_write_png(path.c_str(), img.width, img.height, 1, out_pixels.data(), img.width);
    if (!success)
    {
        throw std::runtime_error("Erro ao salvar a imagem em: " + path);
    }
}

Image convolve3x3(const Image &src, const double k[3][3])
{
    Image dst = src;

    for (int i = 0; i < src.height; ++i)
    {
        for (int j = 0; j < src.width; ++j)
        {
            double acc = 0.0;

            for (int di = -1; di <= 1; ++di)
            {
                for (int dj = -1; dj <= 1; ++dj)
                {
                    int ii = std::clamp(i + di, 0, src.height - 1);
                    int jj = std::clamp(j + dj, 0, src.width - 1);
                    acc += src.at(ii, jj) * k[di + 1][dj + 1];
                }
            }

            dst.at(i, j) = acc;
        }
    }

    return dst;
}

void gaussianKernel(double k[3][3])
{
    double base[3][3] = {{1, 2, 1}, {2, 4, 2}, {1, 2, 1}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            k[i][j] = base[i][j] / 16.0;
}

void sobelX(double k[3][3])
{
    double base[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            k[i][j] = base[i][j];
}

void sobelY(double k[3][3])
{
    double base[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            k[i][j] = base[i][j];
}

void laplaceKernel(double k[3][3])
{
    double base[3][3] = {{0, 1, 0}, {1, -4, 1}, {0, 1, 0}};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            k[i][j] = base[i][j];
}

Image alg1_sobel(const Image &original, double threshold)
{
    double gK[3][3], sx[3][3], sy[3][3];
    gaussianKernel(gK);
    sobelX(sx);
    sobelY(sy);

    Image smooth = convolve3x3(original, gK);

    Image A = convolve3x3(smooth, sx);

    Image B = convolve3x3(smooth, sy);

    Image C = original;
    C.pixels.assign(static_cast<size_t>(original.width) * original.height, 0.0);

    for (int i = 0; i < original.height; ++i)
    {
        for (int j = 0; j < original.width; ++j)
        {
            double ax = A.at(i, j);
            double by = B.at(i, j);
            C.at(i, j) = std::sqrt(ax * ax + by * by);
        }
    }

    Image D = original;
    D.pixels.assign(static_cast<size_t>(original.width) * original.height, 0.0);

    for (int i = 0; i < original.height; ++i)
    {
        for (int j = 0; j < original.width; ++j)
        {
            D.at(i, j) = (C.at(i, j) > threshold) ? 1.0 : 0.0;
        }
    }

    return D;
}

Image alg2_laplace(const Image &original, double tolerance)
{
    double gK[3][3], lK[3][3];
    gaussianKernel(gK);
    laplaceKernel(lK);

    Image smooth = convolve3x3(original, gK);

    Image A = convolve3x3(smooth, lK);

    Image B = original;
    B.pixels.assign(static_cast<size_t>(original.width) * original.height, 0.0);

    for (int i = 0; i < original.height; ++i)
    {
        for (int j = 0; j < original.width; ++j)
        {
            if (std::abs(A.at(i, j)) <= tolerance)
                B.at(i, j) = 0.0;
            else
                B.at(i, j) = 1.0;
        }
    }

    return B;
}

int main(int argc, char **argv)
{
    if (argc < 4)
    {
        std::cerr << "Uso: " << argv[0]
                  << " entrada.png saida_sobel.png saida_laplace.png [threshold=35.0] [tolerancia=0.0001]\n";
        return 1;
    }

    try
    {
        Image img = readRealImage(argv[1]);

        double threshold = (argc > 4) ? std::stod(argv[4]) : 35.0;
        double tolerance = (argc > 5) ? std::stod(argv[5]) : 3.0;

        Image sobel = alg1_sobel(img, threshold);
        Image laplace = alg2_laplace(img, tolerance);

        writeRealImage(argv[2], sobel, true);
        writeRealImage(argv[3], laplace, true);

        std::cout << "Processamento concluído com sucesso.\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "[ERRO] " << e.what() << "\n";
        return 1;
    }

    return 0;
}