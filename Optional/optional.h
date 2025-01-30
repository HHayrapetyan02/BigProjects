#ifndef OPTIONAL_H
#define OPTIONAL_H

#include <iostream>

class BadOptionalAccess : public std::exception {
 public:
  const char* what() const noexcept override {
    return "Bad Optional Access";
  }
};

template <class T>
class Optional {
  bool has_value_ = false;
  char memory_[sizeof(T)];

 public:
  Optional();
  Optional(const Optional<T>& optional);
  Optional(Optional<T>&&) noexcept;
  Optional(const T&);  // NOLINT
  Optional(T&&);       // NOLINT
  ~Optional();
  Optional<T>& operator=(Optional<T>&&) noexcept;
  Optional<T>& operator=(const Optional<T>&);
  Optional<T>& operator=(const T&);
  Optional<T>& operator=(T&&);
  bool HasValue() const noexcept;
  explicit operator bool() const noexcept;
  T& Value();
  const T& Value() const;
  const T& operator*() const;
  T& operator*();
  template <class... Args>
  T& Emplace(Args&&...);
  void Reset();
  void Swap(Optional<T>&);
};

template <class T>
void Optional<T>::Swap(Optional<T>& optional) {
  Optional tmp = std::move(*this);
  *this = std::move(optional);
  optional = std::move(tmp);
}

template <class T>
Optional<T>::Optional() : has_value_(false) {
}

template <class T>
Optional<T>::Optional(const Optional<T>& optional) : has_value_(optional.has_value_) {
  if (has_value_) {
    new (memory_) T(*reinterpret_cast<const T*>(optional.memory_));
  }
}

template <class T>
Optional<T>::Optional(Optional<T>&& optional) noexcept : has_value_(optional.has_value_) {
  if (has_value_) {
    new (memory_) T(std::move(*reinterpret_cast<T*>(optional.memory_)));
  }
}

template <class T>
Optional<T>::Optional(const T& obj) {
  new (memory_) T(obj);
  has_value_ = true;
}

template <class T>
Optional<T>::Optional(T&& obj) {
  new (memory_) T(std::move(obj));
  has_value_ = true;
}

template <class T>
Optional<T>::~Optional() {
  Reset();
}

template <class T>
void Optional<T>::Reset() {
  if (has_value_) {
    has_value_ = false;
    reinterpret_cast<T*>(memory_)->~T();
  }
}

template <class T>
Optional<T>& Optional<T>::operator=(Optional<T>&& optional) noexcept {
  if (this != &optional) {
    if (optional.has_value_) {
      if (has_value_) {
        *reinterpret_cast<T*>(memory_) = std::move(*reinterpret_cast<T*>(optional.memory_));
      } else {
        new (memory_) T(std::move(*reinterpret_cast<T*>(optional.memory_)));
        has_value_ = true;
      }
    } else {
      Reset();
    }
  }
  return *this;
}

template <class T>
Optional<T>& Optional<T>::operator=(const Optional<T>& optional) {
  if (this != &optional) {
    if (optional.has_value_) {
      if (has_value_) {
        *reinterpret_cast<T*>(memory_) = *reinterpret_cast<const T*>(optional.memory_);
      } else {
        new (memory_) T(*reinterpret_cast<const T*>(optional.memory_));
        has_value_ = true;
      }
    } else {
      Reset();
    }
  }
  return *this;
}

template <class T>
bool Optional<T>::HasValue() const noexcept {
  return has_value_;
}

template <class T>
Optional<T>::operator bool() const noexcept {
  return has_value_;
}

template <class T>
Optional<T>& Optional<T>::operator=(const T& val) {
  if (has_value_) {
    *reinterpret_cast<T*>(memory_) = val;
  } else {
    new (memory_) T(val);
    has_value_ = true;
  }
  return *this;
}

template <class T>
Optional<T>& Optional<T>::operator=(T&& val) {
  if (has_value_) {
    *reinterpret_cast<T*>(memory_) = std::move(val);
  } else {
    new (memory_) T(std::move(val));
    has_value_ = true;
  }
  return *this;
}

template <class T>
T& Optional<T>::Value() {
  if (has_value_) {
    return *reinterpret_cast<T*>(memory_);
  }
  throw BadOptionalAccess{};
}

template <class T>
const T& Optional<T>::Value() const {
  if (has_value_) {
    return *reinterpret_cast<const T*>(memory_);
  }
  throw BadOptionalAccess{};
}

template <class T>
T& Optional<T>::operator*() {
  return *reinterpret_cast<T*>(memory_);
}

template <class T>
const T& Optional<T>::operator*() const {
  return *reinterpret_cast<const T*>(memory_);
}

template <class T>
template <class... Args>
T& Optional<T>::Emplace(Args&&... args) {
  Reset();
  has_value_ = true;
  return *(new (memory_) T(std::forward<Args>(args)...));
}
#endif
