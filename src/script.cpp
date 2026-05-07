#include "script.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <stdexcept>

namespace {

struct InitBlock {
    int xb = 0;
    int xe = 0;
    int yb = 0;
    int ye = 0;
    double rho = 0.0;
    double u_vel = 0.0;
    double v_vel = 0.0;
    double p = 0.0;
    double w = 0.0;
    bool has_w = false;
};

std::vector<InitBlock> parse_blocks(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open " + filename);
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();

    std::regex block_pattern(R"(BLOCK_\d+([\s\S]*?)(?=BLOCK_\d+|$))");
    std::regex kv_pattern(R"(([\w_]+)\s*[:=]\s*([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?))");

    auto blocks_begin = std::sregex_iterator(content.begin(), content.end(), block_pattern);
    auto blocks_end = std::sregex_iterator();
    std::vector<InitBlock> blocks;

    for (auto i = blocks_begin; i != blocks_end; ++i) {
        std::string block_body = (*i)[1].str();
        InitBlock block;

        auto kv_begin = std::sregex_iterator(block_body.begin(), block_body.end(), kv_pattern);
        for (auto j = kv_begin; j != std::sregex_iterator(); ++j) {
            std::string key = (*j)[1].str();
            double val = std::stod((*j)[2].str());
            if (key == "x_begin") block.xb = static_cast<int>(val);
            else if (key == "x_end") block.xe = static_cast<int>(val);
            else if (key == "y_begin") block.yb = static_cast<int>(val);
            else if (key == "y_end") block.ye = static_cast<int>(val);
            else if (key == "density") block.rho = val;
            else if (key == "velocity_u") block.u_vel = val;
            else if (key == "velocity_v") block.v_vel = val;
            else if (key == "pressure") block.p = val;
            else if (key == "mass_fraction_w" || key == "w") {
                block.w = val;
                block.has_w = true;
            }
        }
        blocks.push_back(block);
    }

    return blocks;
}

void fill_block(std::vector<std::vector<std::vector<double>>>& u,
                const InitBlock& block,
                int offset_x,
                int offset_y,
                double gamma) {
    const double rhou = block.rho * block.u_vel;
    const double rhov = block.rho * block.v_vel;
    const double energy = block.p / (gamma - 1.0)
        + 0.5 * block.rho * (block.u_vel * block.u_vel + block.v_vel * block.v_vel);

    for (int x = block.xb; x <= block.xe; ++x) {
        for (int y = block.yb; y <= block.ye; ++y) {
            int tx = x + offset_x;
            int ty = y + offset_y;
            if (tx >= 0 && ty >= 0
                && tx < static_cast<int>(u.size())
                && ty < static_cast<int>(u[0].size())) {
                u[tx][ty][0] = block.rho;
                u[tx][ty][1] = rhou;
                u[tx][ty][2] = rhov;
                u[tx][ty][3] = energy;
            }
        }
    }
}

} // namespace

void parseAndInitialize(std::vector<std::vector<std::vector<double>>>& u,
                        const std::string& filename,
                        int offset_x,
                        int offset_y,
                        double gamma) {
    const std::vector<InitBlock> blocks = parse_blocks(filename);
    for (const InitBlock& block : blocks) {
        fill_block(u, block, offset_x, offset_y, gamma);
    }
}

void parseAndInitializeMader(std::vector<std::vector<std::vector<double>>>& u,
                             ScalarField& w,
                             const std::string& filename,
                             int offset_x,
                             int offset_y,
                             double gamma,
                             double default_w) {
    const std::vector<InitBlock> blocks = parse_blocks(filename);
    for (const InitBlock& block : blocks) {
        fill_block(u, block, offset_x, offset_y, gamma);
        const double block_w = block.has_w ? block.w : default_w;
        for (int x = block.xb; x <= block.xe; ++x) {
            for (int y = block.yb; y <= block.ye; ++y) {
                int tx = x + offset_x;
                int ty = y + offset_y;
                if (tx >= 0 && ty >= 0
                    && tx < static_cast<int>(w.size())
                    && ty < static_cast<int>(w[0].size())) {
                    w[tx][ty] = block_w;
                }
            }
        }
    }
}
