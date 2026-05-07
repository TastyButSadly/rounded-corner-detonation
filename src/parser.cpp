#include "parser.h"


SectionedIniParser::SectionedIniParser(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    config.clear();
    std::string current_section = "global"; // секция по умолчанию
    
    std::string line;
    while (std::getline(file, line)) {
        // Безопасный тримминг
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        size_t end = line.find_last_not_of(" \t\r\n");
        std::string trimmed = line.substr(start, end - start + 1);
        
        if (trimmed.empty()) continue;
        
        // Комментарии
        if (trimmed[0] == ';' || trimmed[0] == '#') continue;
        
        // Секции
        if (trimmed[0] == '[') {
            size_t end_bracket = trimmed.find(']');
            if (end_bracket != std::string::npos && end_bracket > 1) {
                current_section = trimmed.substr(1, end_bracket - 1);
                // Дополнительный тримминг
                current_section.erase(0, current_section.find_first_not_of(" \t"));
                current_section.erase(current_section.find_last_not_of(" \t") + 1);
                if (!current_section.empty()) {
                    config[current_section] = std::map<std::string, std::string>();
                }
            }
            continue;
        }
        
        // Ключ-значение
        size_t equals_pos = trimmed.find('=');
        if (equals_pos != std::string::npos && equals_pos > 0) {
            std::string key = trimmed.substr(0, equals_pos);
            std::string value = trimmed.substr(equals_pos + 1);
            
            // Тримминг ключа и значения
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);
            
            if (!key.empty() && !current_section.empty()) {
                config[current_section][key] = value;
            }
        }
    }
}


std::string SectionedIniParser::getString(const std::string& section, const std::string& key) const {
    // Тримминг входных параметров
    std::string section_trimmed = trim(section);
    std::string key_trimmed = trim(key);
    
    // Поиск секции
    auto sec_it = config.find(section_trimmed);
    if (sec_it == config.end()) {
        throw std::runtime_error("Section '" + section_trimmed + "' not found in ini file.");
    }
    
    // Поиск ключа в секции
    auto key_it = sec_it->second.find(key_trimmed);
    if (key_it == sec_it->second.end()) {
        throw std::runtime_error("Key '" + key_trimmed + "' not found in section '" + section_trimmed + "'.");
    }
    // if (key_it == sec_it->second.end()) {
    //     // Отладочный вывод доступных ключей
    //     std::string available_keys;
    //     for (const auto& k : sec_it->second) {
    //         available_keys += "'" + k.first + "' ";
    //     }
    //     throw std::runtime_error("Key '" + key_trimmed + "' not found in section '" + section_trimmed + "'. Available keys: " + available_keys);
    // }
    
    return key_it->second;
}

   
std::string SectionedIniParser::trim(const std::string& s) const{
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}


IniParser::IniParser(const std::string& filename) {
    std::ifstream fin(filename);
    if (!fin.is_open())
        throw std::runtime_error("Cannot open ini file: " + filename);

    std::string line;
    while (std::getline(fin, line)) {
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
            line = line.substr(0, comment_pos);

        trim(line);
        if (line.empty()) continue;

        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        std::string key = line.substr(0, eq);
        std::string value = line.substr(eq + 1);
        trim(key);
        trim(value);
        params[key] = value;
    }
}

std::string IniParser::getString(const std::string& key) const {
    auto it = params.find(key);
    if (it == params.end())
        throw std::runtime_error("Key '" + key + "' not found.");
    return it->second;
}


void IniParser::trim(std::string& s) {
    size_t first = s.find_first_not_of(" \t");
    size_t last = s.find_last_not_of(" \t");
    if (first == std::string::npos || last == std::string::npos) {
        s.clear();
    } else {
        s = s.substr(first, last - first + 1);
    }
}

