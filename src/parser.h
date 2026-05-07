#ifndef PARSER
#define PARSER

#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

class SectionedIniParser {
public:
    explicit SectionedIniParser(const std::string& filename);

    std::string getString(const std::string& section, const std::string& key) const;

    double getDouble(const std::string& section, const std::string& key) const {
        return std::stod(getString(section, key));
    }

    int getInt(const std::string& section, const std::string& key) const {
        return std::stoi(getString(section, key));
    }

    bool hasSection(const std::string& section) const {
        return config.count(section) > 0;
    }

private:
    std::map<std::string, std::map<std::string, std::string> > config;

    std::string trim(const std::string& s) const;
};


class IniParser {
public:
    explicit IniParser(const std::string& filename);

    std::string getString(const std::string& key) const;

    double getDouble(const std::string& key) const {
        return std::stod(getString(key));
    }

    int getInt(const std::string& key) const {
        return std::stoi(getString(key));
    }

private:
    std::map<std::string, std::string> params;

    static void trim(std::string& s);
};


#endif
