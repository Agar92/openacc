#pragma once
#ifndef T3PROCESSIMPLFROMFSANDCS_H
#define T3PROCESSIMPLFROMFSANDCS_H

#include <memory>

namespace t3 {

template <typename CSImpl, typename FSImpl> class ProcessImplFromCSandFS;

template <typename CSImpl, typename FSImpl> class ProcessImplMakerMedium;

template <typename CSImpl, typename FSImpl, typename... Args>
ProcessImplMakerMedium<CSImpl, FSImpl>
makeProcessImplFromCSandFS(Args &&... args);

template <typename CSImpl, typename FSImpl> class ProcessImplFromCSandFS {
public:
  using CS_t = CSImpl;
  using FS_t = FSImpl;
  ProcessImplFromCSandFS(CSImpl &&cs, FSImpl &&fs)
      : fCS(std::make_unique<CSImpl>(std::move(cs))),
        fFS(std::make_unique<FSImpl>(std::move(fs))) {}
  template <typename Floating, typename... Args>
  auto GetCS(Floating e, Args &&... args) const;
  template <typename... Args> auto GetFS(Args &&... args) const;

protected:
  std::unique_ptr<CSImpl> fCS;
  std::unique_ptr<FSImpl> fFS;
};

template <typename CSImpl, typename FSImpl>
template <typename Floating, typename... Args>
auto ProcessImplFromCSandFS<CSImpl, FSImpl>::GetCS(Floating e,
                                                   Args &&... args) const {
  return fCS->GetCS(e, std::forward<Args>(args)...);
}

template <typename CSImpl, typename FSImpl>
template <typename... Args>
auto ProcessImplFromCSandFS<CSImpl, FSImpl>::GetFS(Args &&... args) const {
  return fFS->GetFS(std::forward<Args>(args)...);
}

template <typename CSImpl, typename FSImpl> class ProcessImplMakerMedium {
public:
  template <typename... Args>
  ProcessImplMakerMedium(Args &&... args) : fCS(std::forward<Args>(args)...) {}
  template <typename... Args>
  ProcessImplFromCSandFS<CSImpl, FSImpl> operator()(Args &&... args) {
    return ProcessImplFromCSandFS<CSImpl, FSImpl>(
        fCS, FSImpl(std::forward<Args>(args)...));
  }

private:
  CSImpl fCS;
};

template <typename CSImpl, typename FSImpl, typename... Args>
ProcessImplMakerMedium<CSImpl, FSImpl>
makeProcessImplFromCSandFS(Args &&... args) {
  return ProcessImplMakerMedium<CSImpl, FSImpl>(std::forward<Args>(args)...);
}

} // namespace t3

#endif // T3PROCESSIMPLFROMFSANDCS_H
