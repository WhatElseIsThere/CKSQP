#include <algorithm>
#include <vector>

namespace CKSQP {

using VecD = std::vector<double>;

class Func {
public:
  Func(double alpha, double lambda) : alpha_(alpha), lambda_(lambda) {}

  double Alpha() const {
    return alpha_;
  }

  double Lambda() const {
    return lambda_;
  }

  double Eval(double lambda) const {
    return alpha_ * std::max(lambda - lambda_, .0);
  }

  bool operator<(const Func &other) const {
    return lambda_ < other.Lambda();
  }

private:
  double alpha_, lambda_;
};


class FuncList {
public:
  void AddFunc(Func func) { list_.push_back(func); }

  double Eval(double lambda, size_t l, size_t r) const {
    double res = 0;
    for (auto i = l; i < r; ++i) {
      auto func = list_[i];
      res += func.Eval(lambda);
    }
    return res;
  }

  double Eval(double lambda) const { return Eval(lambda, 0, list_.size()); }

  const std::vector<Func> &Vec() const { return list_; }

  std::vector<Func> &Vec() { return list_; }

private:
  std::vector<Func> list_;
};


struct Equation {
  FuncList list;
  double c;
};


class Problem {
public:
  Problem(size_t n, VecD h, VecD l, VecD r, double b)
      : n_(n), h_(h), l_(l), r_(r), b_(b) {}

  size_t N() const {
    return n_;
  }

  Equation BuildEquation() const {
    FuncList list;
    double c = b_;
    for (int i = 0; i < n_; ++i) {
      list.AddFunc(Func(+1 / h_[i], h_[i] * l_[i]));
      list.AddFunc(Func(-1 / h_[i], h_[i] * r_[i]));
      c -= l_[i];
    }
    return {list, c};
  }

  // std::vector<double> BuildSolution(double lambda) const {
  //   std::vector<double> res;
  //   for (size_t i = 0; i < h_.size(); ++i) {
  //     double x = l_[i] + (std::max(lambda - h_[i] * l_[i], 0.0) -
  //                         std::max(lambda - h_[i] * r_[i], 0.0)) /
  //                            h_[i];
  //     res.push_back(x);
  //   }
  //   return res;
  // }

private:
  size_t n_;
  std::vector<double> h_, l_, r_;
  double b_;
};


class Solver {
public:
  virtual double solve(const Problem &P) = 0;
};


class LinSolver : Solver {
public:
  double solve(const Problem &P) {
    auto n = P.N();

    auto eq = P.BuildEquation();
    auto &list = eq.list;
    auto c = eq.c;

    size_t lo = 0, hi = 2 * n;
    double a = 0, b = 0;

    while (hi - lo > 1) {
      size_t mid = lo + (hi - lo) / 2;

      auto &vec = list.Vec();
      std::nth_element(vec.begin() + lo, vec.begin() + mid, vec.begin() + hi);

      auto lambda = vec[mid].Lambda();
      auto value = a * lambda + b + list.Eval(lambda, lo, hi);

      if (value <= c) {
        for (auto i = lo; i < mid; ++i) {
          auto func = vec[i];
          a += func.Alpha();
          b -= func.Alpha() * func.Lambda();
        }
        lo = mid;
      } else {
        hi = mid;
      }
    }

    return list.Vec()[lo].Lambda();
  }
};


class BinSolver : Solver {
public:
  double solve(const Problem &P) {
    auto n = P.N();

    auto eq = P.BuildEquation();
    auto &list = eq.list;
    auto c = eq.c;

    size_t lo = 0, hi = 2 * n;
    while (hi - lo > 1) {
      size_t mid = lo + (hi - lo) / 2;
      auto &vec = list.Vec();
      std::nth_element(vec.begin() + lo, vec.begin() + mid, vec.begin() + hi);

      auto lambda = vec[mid].Lambda();
      auto value = list.Eval(lambda);

      if (value <= c) {
        lo = mid;
      } else {
        hi = mid;
      }
    }

    return list.Vec()[lo].Lambda();
  }
};


class DumbSolver : Solver {
public:
  double solve(const Problem &P) {
    auto n = P.N();

    auto eq = P.BuildEquation();
    const auto &list = eq.list;
    const auto &vec = list.Vec();
    auto c = eq.c;

    double ans = std::numeric_limits<double>::lowest();
    for (auto func : vec) {
      auto lambda = func.Lambda();
      if (list.Eval(lambda) <= c) {
        ans = std::max(ans, lambda);
      }
    }

    return ans;
  }
};

} // namespace CKSQP
