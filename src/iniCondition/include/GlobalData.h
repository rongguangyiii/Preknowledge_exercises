/*---------------------------------------------------------------------------
	Pre knowledge exercises
	Master the basic skills of equation writing and calculation
	Copyright (C) ,Since 2024
-------------------------------------------------------------------------------
License
!	@file		FLOWPARAMETER.h
!	@brief	initial flow information.
!	@author	liu guanying.
\*---------------------------------------------------------------------------*/
#ifndef FLOWPARAMETER_H
#define FLOWPARAMETER_H
#include <unordered_map>
#include <variant>
#include <string>

class GlobalData
{
public:
    using dataVariant = std::variant<std::string, int, double>;
    GlobalData(const GlobalData&) = delete;
    GlobalData& operator=(const GlobalData&) = delete;
    static GlobalData& Init();
    static int GetInt(const std::string& dataName);
    static double GetDouble(const std::string& dataName);
    static std::string GetString(const std::string& dataName);
    static void Update(const std::string& dataName, const dataVariant& dataValue);

private:
    GlobalData() {}
    ~GlobalData() {}
    static const dataVariant& Get(const std::string& dataName);
    static std::unordered_map<std::string, dataVariant> dataMap_;
};


#endif