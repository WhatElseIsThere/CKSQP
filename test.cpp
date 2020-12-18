#include <bits/stdc++.h>
#include "cksqp.hpp"

using namespace std;
using namespace CKSQP;

mt19937 rng(random_device{}());

Problem GenProblem(size_t n) {
  uniform_real_distribution<> h_rd(0.1, 1.0);
  vector<double> h(n);
  for (int i = 0; i < n; ++i) {
    h[i] = h_rd(rng);
  }

  uniform_real_distribution<> lr_rd(-1.0, +1.0);
  vector<double> l(n), r(n);
  for (int i = 0; i < n; ++i) {
    l[i] = lr_rd(rng);
    r[i] = lr_rd(rng);
    if (l[i] > r[i]) swap(l[i], r[i]);
  }

  auto l_sum = accumulate(l.begin(), l.end(), .0);
  auto r_sum = accumulate(r.begin(), r.end(), .0);

  uniform_real_distribution<> b_rd(l_sum, r_sum);
  double b = b_rd(rng);

  return Problem(n, h, l, r, b);
}

int main() {
  LinSolver ls;
  BinSolver bs;
  DumbSolver ds;

  size_t n;
  cin >> n;

  const size_t its = 100;

  double lin_time = .0, bin_time = .0, dumb_time = .0;

  for (int it = 0; it < its; ++it) {
    auto P = GenProblem(n);
    vector<double> answers;

    clock_t start, finish;

    start = clock();
    answers.push_back(ls.solve(P));
    finish = clock();
    lin_time += double(finish - start) / CLOCKS_PER_SEC;

    start = clock();
    answers.push_back(bs.solve(P));
    finish = clock();
    bin_time += double(finish - start) / CLOCKS_PER_SEC;

    // start = clock();
    // answers.push_back(ds.solve(P));
    // finish = clock();
    // dumb_time += double(finish - start) / CLOCKS_PER_SEC;

    auto min_ans = *min_element(answers.begin(), answers.end());
    auto max_ans = *max_element(answers.begin(), answers.end());

    if (min_ans == max_ans) {
      cout << it + 1 << " : \033[32m[OK]\033[0m" << endl;
    } else {
      cout << it + 1 << " : \033[31m[ERROR]\033[0m" << endl;
      for (auto ans : answers) { cout << ans << ' '; }
      cout << endl;
    }
  }

  cout << setprecision(10) << fixed;
  cout << "BinSolver time: " << ceil(bin_time / its * 1000) << endl;
  cout << "LinSolver time: " << ceil(lin_time / its * 1000) << endl;
  // cout << "DumbSolver time: " << dumb_time / its << endl;
}
