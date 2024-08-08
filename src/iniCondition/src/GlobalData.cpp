#include "iniCondition/include/GlobalData.h"
#include "utility/include/log.h"
#include <fstream>

GlobalData& GlobalData::Init(const std::string& fileName)
{
	static GlobalData data;
	data.readconfig(fileName);
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
}

void GlobalData::readconfig(const std::string& fileName)
{
	std::ifstream fin("../../config/" + fileName);
	std::string line;
	std::string dataType;
	std::string dataName;
	std::string dataValue;

	while (std::getline(fin, line))
	{
		std::string separator = " =\r\n\t#$,;\"";
		if (line.empty())
			continue;
		size_t id = line.find_first_not_of(separator);
		if (id == std::string::npos)
			continue;
		//删除前方没用的信息
		line.erase(line.begin(), line.begin() + id);
		//注释行
		if (line[0] == '!' || line[0] == '！' || line[0] == '#' || line[0] == '/')
			continue;
		id = line.find_first_of(separator);
		dataType = line.substr(0, id);
		line.erase(0, id);
		line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
		if (line.empty())
			continue;
		id = line.find_first_of("=");
		if (id == std::string::npos)
			continue;
		dataName = line.substr(0, id);
		dataValue = line.substr(id + 1);
		if (dataType == "string")
			GlobalData::Update(dataName, dataValue);
		else if (dataType == "double")
			GlobalData::Update(dataName, stod(dataValue));
		else if (dataType == "int")
			GlobalData::Update(dataName, stoi(dataValue));
		else
		{
			spdlog::warn("Unspported Data Type:{}, Name:{}, Value:{}", dataType, dataName, dataValue);
		}
	}
}
