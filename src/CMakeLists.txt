set(CMAKE_CXX_STANDARD 17)

#自动查找指定目录下的所有源文件
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/Coordinate/src SRC_FILE)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/flowControl/src SRC_FILE)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/mesh/src SRC_FILE)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/iniCondition/src SRC_FILE)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/utility/src SRC_FILE)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/solver/src SRC_FILE)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/outfile/src SRC_FILE)

# 将子目录中的源文件添加到父项目的源文件列表中
# 使用PARENT_SCOPE选项与set命令一起时，实际上是在设置一个上下文都可见的变量。
set(SRC_FILE ${SRC_FILE} PARENT_SCOPE)

