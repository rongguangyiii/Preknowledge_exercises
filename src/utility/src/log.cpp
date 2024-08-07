#include "utility/include/log.h"
#include <cstdio>
#include <chrono>
Log& Log::Satrt()
{
	static Log log;
	return log;
}
Log::Log()
{
	auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
	console_sink->set_level(spdlog::level::debug);
	console_sink->set_pattern("[%D %H:%M:%S] [%^%l%$] %v");
	auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>("log/log.txt", true);
	file_sink->set_level(spdlog::level::trace);
	file_sink->set_pattern("[%D %H:%M:%S] [%^%l%$] %v");
	spdlog::logger logger("multi_sink", { console_sink, file_sink });
	logger.set_level(spdlog::level::trace);
	logger.info("Start Logger, Log file path: ./{}", file_sink->filename());
	spdlog::set_default_logger(std::make_shared<spdlog::logger>(logger));
	spdlog::flush_every(std::chrono::seconds(15));//ÿ��15sˢ��һ��log�ļ�
}

Log::~Log()
{
	spdlog::drop_all();
	spdlog::shutdown();
}