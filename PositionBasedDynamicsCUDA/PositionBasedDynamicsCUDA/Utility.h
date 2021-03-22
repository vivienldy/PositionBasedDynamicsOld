#pragma once

#include <vector>
#include <string>

namespace UT
{
	inline void splitString(const std::string& src, std::vector<std::string>& v, const std::string& split)
	{
		std::string::size_type pos1, pos2;
		pos2 = src.find(split);
		pos1 = 0;
		while (std::string::npos != pos2)
		{
			v.push_back(src.substr(pos1, pos2 - pos1));
			pos1 = pos2 + split.size();
			pos2 = src.find(split, pos1);
		}
		if (pos1 != src.length())
			v.push_back(src.substr(pos1));
	}

	inline std::vector<std::string> splitString(const std::string& src, const std::string& split)
	{
		std::vector<std::string> _ret = std::vector<std::string>();
		splitString(src, _ret, split);
		return _ret;
	}
}
