#ifndef SCRIPT
#define SCRIPT

#include <vector>
#include <string>

#include "mader.h"

// Добавляем параметры для учета фиктивных ячеек и физических констант
void parseAndInitialize(std::vector<std::vector<std::vector<double>>>& u, 
                        const std::string& filename, 
                        int offset_x, 
                        int offset_y, 
                        double gamma);

void parseAndInitializeMader(std::vector<std::vector<std::vector<double>>>& u,
                             ScalarField& w,
                             const std::string& filename,
                             int offset_x,
                             int offset_y,
                             double gamma,
                             double default_w);

#endif
