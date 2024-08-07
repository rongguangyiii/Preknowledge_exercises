#include "iniCondition/include/GlobalData.h"
#include "utility/include/log.h"
#include <fstream>

GlobalData& GlobalData::Init()
{
	static GlobalData data;
    return data;
}

std::unordered_map<std::string, GlobalData::dataVariant> GlobalData::dataMap_;

const GlobalData::dataVariant& GlobalData::Get(const std::string& dataName)
{
    auto it = dataMap_.find(dataName);
    if (it != dataMap_.end())
        return it->second;
    else
    {
        spdlog::warn("Not Found: {}, in Global Data!", dataName);
        spdlog::warn("Exit Now!");
        exit(0);
    }
}

int GlobalData::GetInt(const std::string& dataName)
{
    const dataVariant& value = Get(dataName);
    if (std::holds_alternative<int>(value))
        return std::get<int>(value);
    else
        throw std::runtime_error("Data type is not int for key: " + dataName);
}

double GlobalData::GetDouble(const std::string& dataName)
{
    const dataVariant& value = Get(dataName);
    //spdlog::info("Retrieved {} with type index: {}", dataName, value.index()); // debug info
    if (std::holds_alternative<double>(value))
        return std::get<double>(value);
    else
        throw std::runtime_error("Data type is not double for key: " + dataName);
}

std::string GlobalData::GetString(const std::string& dataName)
{
    const dataVariant& value = Get(dataName);
    if (std::holds_alternative<std::string>(value))
        return std::get<std::string>(value);
    else
        throw std::runtime_error("Data type is not string for key: " + dataName);
}

void GlobalData::Update(const std::string& varName, const dataVariant& varValue)
{
    dataMap_[varName] = varValue;
    //spdlog::info("Updated {} with value: {}", varName, varValue.index()); // debug info
}
