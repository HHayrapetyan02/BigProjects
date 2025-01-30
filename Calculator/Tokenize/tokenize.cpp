#include "tokenize.h"

bool PlusToken::operator==(const PlusToken&) const {
  return true;
}

bool MinusToken::operator==(const MinusToken&) const {
  return true;
}

bool MultiplyToken::operator==(const MultiplyToken&) const {
  return true;
}

bool DivideToken::operator==(const DivideToken&) const {
  return true;
}

bool ResidualToken::operator==(const ResidualToken&) const {
  return true;
}

bool OpeningBracketToken::operator==(const OpeningBracketToken&) const {
  return true;
}

bool ClosingBracketToken::operator==(const ClosingBracketToken&) const {
  return true;
}

bool SqrToken::operator==(const SqrToken&) const {
  return true;
}

bool MaxToken::operator==(const MaxToken&) const {
  return true;
}

bool MinToken::operator==(const MinToken&) const {
  return true;
}

bool AbsToken::operator==(const AbsToken&) const {
  return true;
}

UnknownToken::UnknownToken(std::string str) : value(str) {
}

bool UnknownToken::operator==(const UnknownToken& token) const {
  return value == token.value;
}

NumberToken::NumberToken(int value) : value(value) {
}

bool NumberToken::operator==(const NumberToken& token) const {
  return value == token.value;
}

std::vector<Token> Tokenize(const std::string_view& str) {
  std::string string = std::string();
  std::vector<std::string> vector_strings;
  for (size_t i = 0; i < str.size(); ++i) {
    if (str[i] == ' ' && string.empty()) {
      continue;
    }
    if (str[i] == ' ' && !string.empty()) {
      vector_strings.push_back(string);
      string = std::string();
    } else {
      string.push_back(str[i]);
      if (i == str.size() - 1) {
        vector_strings.push_back(string);
        string = std::string();
      }
    }
  }
  std::vector<Token> ans;
  Token token;
  for (const auto& now_str : vector_strings) {
    if (now_str == "+") {
      token = PlusToken();
      ans.push_back(token);
    } else if (now_str == "-") {
      token = MinusToken();
      ans.push_back(token);
    } else if (now_str == "*") {
      token = MultiplyToken();
      ans.push_back(token);
    } else if (now_str == "/") {
      token = DivideToken();
      ans.push_back(token);
    } else if (now_str == "%") {
      token = ResidualToken();
      ans.push_back(token);
    } else if (now_str == "(") {
      token = OpeningBracketToken();
      ans.push_back(token);
    } else if (now_str == ")") {
      token = ClosingBracketToken();
      ans.push_back(token);
    } else if (now_str == "sqr") {
      token = SqrToken();
      ans.push_back(token);
    } else if (now_str == "max") {
      token = MaxToken();
      ans.push_back(token);
    } else if (now_str == "min") {
      token = MinToken();
      ans.push_back(token);
    } else if (now_str == "abs") {
      token = AbsToken();
      ans.push_back(token);
    } else {
      bool flag = true;
      for (size_t i = 0; i < now_str.size(); ++i) {
        if (i == 0) {
          flag = now_str[i] == '-' || now_str[i] == '+' || (now_str[i] >= '0' && now_str[i] <= '9');
          if (!flag) {
            break;
          }
        } else {
          flag = now_str[i] >= '0' && now_str[i] <= '9';
          if (!flag) {
            break;
          }
        }
      }
      if (flag) {
        token = NumberToken(std::stoi(now_str));
        ans.push_back(token);
      } else {
        token = UnknownToken(now_str);
        ans.push_back(token);
      }
    }
  }
  return ans;
}
