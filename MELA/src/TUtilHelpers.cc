#include <regex>
#include "TUtilHelpers.hh"


void TUtilHelpers::ExpandEnvironmentVariables(std::string& str){
  static std::regex env("\\$\\{([^}]+)\\}");
  std::smatch match;
  while (std::regex_search(str, match, env)){
    const char* s = getenv(match[1].str().c_str());
    const std::string var(s == NULL ? "" : s);
    str.replace(match[0].first, match[0].second, var);
  }
}
