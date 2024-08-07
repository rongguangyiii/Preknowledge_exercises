/*---------------------------------------------------------------------------
    Pre knowledge exercises
    Master the basic skills of equation writing and calculation
    Copyright (C) ,Since 2024 ,from liuguangying. 
-------------------------------------------------------------------------------*/

#include "flowControl/include/flowControl.h"
#include <iostream>

int main() 
{
    Log::Satrt();
    FlowControl flowControl("config.ini");
    flowControl.preProcess();
    flowControl.start();
    flowControl.postProcess();
    return 0;
}

