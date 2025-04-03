/*
Copyright (C) 2016-2023 Deep Genomics Inc. All Rights Reserved.
*/
#include "gk_assert.h"
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>   // for IsDebuggerPresent
#else
#include <cstring>
#include <fcntl.h>
#include <iterator>
#if defined (__APPLE__)
#include <sys/sysctl.h>
#endif
#include <sys/types.h>
#include <unistd.h>
#endif

BEGIN_NAMESPACE_GK

const char* runtime_error::what() const noexcept
{
	if (!std::empty(buf))
		return buf.c_str();

	const char* const msg = std::runtime_error::what();
	try {
		ccast<decltype(buf)&>(buf) = std::format("{}:{}: {}", file, line, msg);
		return buf.c_str();
	} catch (...) {
		return msg;
	}
}

bool is_debugger_running()
{
#if defined(_WIN32)
	return IsDebuggerPresent() != 0;
#elif defined(__APPLE__)
	// Initialize mib, which tells sysctl the info we want, in this case
	// we're looking for information about a specific process ID.
	int mib[] = { CTL_KERN, KERN_PROC, KERN_PROC_PID, getpid() };
	struct kinfo_proc info;
	size_t size         = sizeof(info);
	info.kp_proc.p_flag = 0;
	sysctl(mib, std::size(mib), &info, &size, nullptr, 0);

	return (info.kp_proc.p_flag & P_TRACED) != 0;
#elif defined(__linux__)
	using namespace std::literals;

	int status_fd = open("/proc/self/status", O_RDONLY);
	if (status_fd == -1)
		return false;
	char buf[1024];
	ssize_t num_read = read(status_fd, buf, sizeof(buf));
	close(status_fd);
	if (num_read <= 0)
		return false;
	buf[num_read] = 0;
	auto tracer = "TracerPid:\t"sv;
	const char* pid = std::strstr(buf, tracer.data());
	if (pid == nullptr)
		return false;
	// Our pid is 0 without a debugger, assume this for any pid starting with 0.
	pid += std::size(tracer);
	return pid < std::end(buf) && *pid != '0';
#else
	return false;
#endif
}

END_NAMESPACE_GK
