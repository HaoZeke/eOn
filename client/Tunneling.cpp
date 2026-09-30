/*
** This file is part of eOn.
**
** SPDX-License-Identifier: BSD-3-Clause
**
** Copyright (c) 2010--present, eOn Development Team
** All rights reserved.
**
** Repo:
** https://github.com/TheochemUI/eOn
*/
#include "eon/Tunneling.h"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace eonc::tunneling {

double massWeightedDistance(const Matter &a, const Matter &b) {
  if (a.numberOfAtoms() != b.numberOfAtoms()) {
    throw std::invalid_argument("the structures hold different atom counts");
  }
  const AtomMatrix dr = a.pbc(b.getPositions() - a.getPositions());
  const auto mass = a.getMasses();
  double sum = 0.0;
  for (long i = 0; i < a.numberOfAtoms(); ++i) {
    if (mass(i) <= 0.0) {
      throw std::invalid_argument(
          "every atom needs a positive mass for a mass-weighted path");
    }
    sum += mass(i) * dr.row(i).squaredNorm();
  }
  return std::sqrt(sum);
}

std::vector<double>
massWeightedPath(const std::vector<std::shared_ptr<Matter>> &band) {
  std::vector<double> s{0.0};
  s.reserve(band.size());
  for (size_t i = 1; i < band.size(); ++i) {
    s.push_back(s.back() + massWeightedDistance(*band[i - 1], *band[i]));
  }
  return s;
}

Profile::Profile(std::vector<double> s, std::vector<double> v)
    : s_(std::move(s)),
      v_(std::move(v)),
      m_(s_.size(), 0.0) {
  const size_t n = s_.size();
  if (n < 2 || v_.size() != n) {
    throw std::invalid_argument(
        "a profile needs matching s and V with two points or more");
  }
  for (size_t k = 1; k < n; ++k) {
    if (!(s_[k] > s_[k - 1])) {
      throw std::invalid_argument(
          "the path coordinate must increase along the band");
    }
  }
  for (size_t k = 1; k + 1 < n; ++k) {
    const double h0 = s_[k] - s_[k - 1];
    const double h1 = s_[k + 1] - s_[k];
    const double d0 = (v_[k] - v_[k - 1]) / h0;
    const double d1 = (v_[k + 1] - v_[k]) / h1;
    if (d0 * d1 <= 0.0) {
      m_[k] = 0.0;
    } else {
      const double w1 = 2.0 * h1 + h0;
      const double w2 = h1 + 2.0 * h0;
      m_[k] = (w1 + w2) / (w1 / d0 + w2 / d1);
    }
  }
}

double Profile::operator()(double x) const {
  x = std::clamp(x, s_.front(), s_.back());
  auto it = std::upper_bound(s_.begin(), s_.end(), x);
  size_t k = static_cast<size_t>(std::distance(s_.begin(), it));
  k = std::clamp<size_t>(k == 0 ? 0 : k - 1, 0, s_.size() - 2);
  const double h = s_[k + 1] - s_[k];
  const double t = (x - s_[k]) / h;
  const double t2 = t * t;
  const double t3 = t2 * t;
  return (2 * t3 - 3 * t2 + 1) * v_[k] + (t3 - 2 * t2 + t) * h * m_[k] +
         (-2 * t3 + 3 * t2) * v_[k + 1] + (t3 - t2) * h * m_[k + 1];
}

double wellCurvature(const Profile &p, bool leftEnd) {
  const auto &s = p.s();
  const auto &v = p.v();
  const size_t n = s.size();
  const double top = *std::max_element(v.begin(), v.end());
  const double floor = leftEnd ? v.front() : v.back();
  const double half = 0.5 * (top - floor);
  // Sums for the normal equations of y = a x^2 + b x^3.
  double s44 = 0, s45 = 0, s55 = 0, sy2 = 0, sy3 = 0;
  size_t used = 0;
  for (size_t j = 1; j < n; ++j) {
    const size_t i = leftEnd ? j : n - 1 - j;
    const double x = leftEnd ? s[i] - s.front() : s.back() - s[i];
    const double y = v[i] - floor;
    if (y > half) {
      break;
    }
    const double x2 = x * x;
    s44 += x2 * x2;
    s45 += x2 * x2 * x;
    s55 += x2 * x2 * x2;
    sy2 += y * x2;
    sy3 += y * x2 * x;
    ++used;
  }
  if (used == 0) {
    // The next image already stands above half the barrier: the parabola
    // through it is all the band says about this well.
    const size_t i = leftEnd ? 1 : n - 2;
    const double x = leftEnd ? s[i] - s.front() : s.back() - s[i];
    return 2.0 * (v[i] - floor) / (x * x);
  }
  if (used == 1) {
    return 2.0 * sy2 / s44;
  }
  const double det = s44 * s55 - s45 * s45;
  const double a = (sy2 * s55 - sy3 * s45) / det;
  return 2.0 * a;
}

double hbarOmega(double curvature) {
  if (!(curvature > 0.0)) {
    throw std::invalid_argument("a well needs a positive curvature");
  }
  return kHbar * std::sqrt(curvature);
}

double wkbAction(const Profile &p, double energy, int points) {
  const double a = p.s().front();
  const double b = p.s().back();
  const double h = (b - a) / (points - 1);
  double sum = 0.0;
  for (int i = 0; i < points; ++i) {
    const double gap = p(a + i * h) - energy;
    const double f = gap > 0.0 ? std::sqrt(2.0 * gap) : 0.0;
    sum += (i == 0 || i == points - 1) ? 0.5 * f : f;
  }
  return sum * h / kHbar;
}

double Splitting::tlsEnergy() const { return std::hypot(delta, delta0); }

Splitting wkbSplitting(const Profile &p, double hwReactant, double hwProduct) {
  const auto &v = p.v();
  Splitting out;
  const double top = *std::max_element(v.begin(), v.end());
  out.delta = v.back() - v.front();
  out.barrier = top - v.front();
  out.hwReactant = hwReactant;
  out.hwProduct = hwProduct;
  out.referenceEnergy =
      std::max(v.front() + 0.5 * hwReactant, v.back() + 0.5 * hwProduct);
  out.action = wkbAction(p, out.referenceEnergy);
  const double hw = std::sqrt(hwReactant * hwProduct);
  out.delta0 = hw / std::numbers::pi * std::exp(-out.action);
  out.deepWells =
      (top - v.front()) > hwReactant && (top - v.back()) > hwProduct;
  return out;
}

Splitting bandSplitting(const std::vector<std::shared_ptr<Matter>> &band,
                        double referenceEnergy) {
  std::vector<double> v;
  v.reserve(band.size());
  for (const auto &image : band) {
    v.push_back(image->getPotentialEnergy() - referenceEnergy);
  }
  const Profile p(massWeightedPath(band), std::move(v));
  return wkbSplitting(p, hbarOmega(wellCurvature(p, true)),
                      hbarOmega(wellCurvature(p, false)));
}

namespace {

/// MatrixXd is row-major; Eigen's decompositions read a column-major
/// triangle, so every symmetric eigenproblem and LU here runs on a copy.
using ColMajorXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

/// Eigenpairs of a symmetric matrix, ascending, decomposed column-major.
Eigen::SelfAdjointEigenSolver<ColMajorXd> symEigen(const MatrixXd &a,
                                                   bool vectors = true) {
  const ColMajorXd s = 0.5 * (a + a.transpose());
  return Eigen::SelfAdjointEigenSolver<ColMajorXd>(
      s, vectors ? Eigen::ComputeEigenvectors : Eigen::EigenvaluesOnly);
}

} // namespace

namespace {

// J = c T (x) I + blockdiag(A_k), T the tridiagonal (2, -1) spring matrix;
// the diagonal blocks hold 2c already, the off-diagonal blocks are -c I.
// Block LU: D_1 = A_1, D_k = A_k - c^2 D_{k-1}^-1.
class BlockChain {
public:
  BlockChain(double c, const std::vector<MatrixXd> &diag)
      : c_(c) {
    lu_.reserve(diag.size());
    for (size_t k = 0; k < diag.size(); ++k) {
      ColMajorXd d = diag[k];
      if (k > 0) {
        d -= c_ * c_ * lu_.back().inverse();
      }
      lu_.emplace_back(std::move(d));
      const auto &f = lu_.back();
      sign_ *= static_cast<int>(std::lround(f.permutationP().determinant()));
      for (long i = 0; i < d.rows(); ++i) {
        const double u = f.matrixLU()(i, i);
        if (u == 0.0) {
          throw std::runtime_error("instanton: singular chain Hessian block");
        }
        sign_ *= u < 0.0 ? -1 : 1;
        logAbsDet_ += std::log(std::abs(u));
      }
    }
  }
  double logAbsDet() const { return logAbsDet_; }
  int sign() const { return sign_; }
  // x = J^-1 b, b and x stacked by bead.
  std::vector<VectorXd> solve(const std::vector<VectorXd> &b) const {
    const size_t m = b.size();
    std::vector<VectorXd> y(m), x(m);
    y[0] = b[0];
    for (size_t k = 1; k < m; ++k) {
      y[k] = b[k] + c_ * lu_[k - 1].solve(y[k - 1]);
    }
    x[m - 1] = lu_[m - 1].solve(y[m - 1]);
    for (size_t k = m - 1; k-- > 0;) {
      x[k] = lu_[k].solve(y[k] + c_ * x[k + 1]);
    }
    return x;
  }

private:
  double c_;
  std::vector<Eigen::PartialPivLU<ColMajorXd>> lu_;
  double logAbsDet_ = 0.0;
  int sign_ = 1;
};

double dot(const std::vector<VectorXd> &a, const std::vector<VectorXd> &b) {
  double s = 0.0;
  for (size_t k = 0; k < a.size(); ++k) {
    s += a[k].dot(b[k]);
  }
  return s;
}

void scale(std::vector<VectorXd> &a, double f) {
  for (auto &v : a) {
    v *= f;
  }
}

// Point at arc-length fraction f along a polyline.
VectorXd alongPolyline(const std::vector<VectorXd> &pts,
                       const std::vector<double> &cum, double f) {
  const double target = f * cum.back();
  const auto it = std::upper_bound(cum.begin(), cum.end(), target);
  const size_t k = std::clamp<size_t>(
      static_cast<size_t>(std::distance(cum.begin(), it)), 1, pts.size() - 1);
  const double seg = cum[k] - cum[k - 1];
  const double t = seg > 0.0 ? (target - cum[k - 1]) / seg : 0.0;
  return pts[k - 1] + std::clamp(t, 0.0, 1.0) * (pts[k] - pts[k - 1]);
}

struct ActionEval {
  double action = 0.0;
  std::vector<double> v;         // interior beads
  std::vector<VectorXd> grad;    // dS/dq over interior beads
  std::vector<VectorXd> potGrad; // dV/dq over interior beads
};

ActionEval evaluateAction(const std::vector<VectorXd> &interior,
                          const VectorXd &start, const VectorXd &end,
                          double vStart, double vEnd, double dtau,
                          const BatchPotential &potential) {
  ActionEval out;
  potential(interior, out.v, out.potGrad);
  const size_t m = interior.size();
  if (out.v.size() != m || out.potGrad.size() != m) {
    throw std::runtime_error("instanton: potential returned the wrong count");
  }
  auto bead = [&](size_t j) -> const VectorXd & {
    return j == 0 ? start : (j == m + 1 ? end : interior[j - 1]);
  };
  double kinetic = 0.0;
  for (size_t j = 0; j <= m; ++j) {
    kinetic += (bead(j + 1) - bead(j)).squaredNorm();
  }
  double pot = 0.5 * (vStart + vEnd);
  for (double vj : out.v) {
    pot += vj;
  }
  out.action = 0.5 * kinetic / dtau + dtau * pot;
  out.grad.resize(m);
  for (size_t j = 1; j <= m; ++j) {
    out.grad[j - 1] = (2.0 * bead(j) - bead(j - 1) - bead(j + 1)) / dtau +
                      dtau * out.potGrad[j - 1];
  }
  return out;
}

double largestBeadNorm(const std::vector<VectorXd> &g) {
  double m = 0.0;
  for (const auto &v : g) {
    m = std::max(m, v.norm());
  }
  return m;
}

} // namespace

double pathOmega(const MatrixXd &hessStart, const MatrixXd &hessEnd,
                 const VectorXd &start, const VectorXd &end) {
  const VectorXd d = (end - start).normalized();
  const double k = std::max(d.dot(hessStart * d), d.dot(hessEnd * d));
  if (!(k > 0.0)) {
    throw std::invalid_argument(
        "pathOmega: no positive curvature along the path at either minimum");
  }
  return std::sqrt(k);
}

Instanton optimizeInstanton(const VectorXd &start, const VectorXd &end,
                            double betaHbar, std::vector<VectorXd> guess,
                            const BatchPotential &potential,
                            const InstantonOptions &options) {
  const long P = options.beads;
  if (P < 4 || !(betaHbar > 0.0) || start.size() != end.size()) {
    throw std::invalid_argument("optimizeInstanton: need P >= 4, beta hbar > 0 "
                                "and ends of one dimension");
  }
  Instanton inst;
  inst.betaHbar = betaHbar;
  inst.dtau = betaHbar / static_cast<double>(P);
  const double dtau = inst.dtau;

  std::vector<double> vEnds;
  std::vector<VectorXd> gEnds;
  potential({start, end}, vEnds, gEnds);
  if (vEnds.size() != 2) {
    throw std::runtime_error("instanton: potential returned the wrong count");
  }
  inst.asymmetry = vEnds[1] - vEnds[0];

  // Beads along the guess (or the straight line) on a tanh kink centred at
  // beta hbar / 2 whose width follows the harmonic decay of a well.
  if (guess.size() < 2) {
    guess = {start, end};
  }
  std::vector<double> cum(guess.size(), 0.0);
  for (size_t k = 1; k < guess.size(); ++k) {
    cum[k] = cum[k - 1] + (guess[k] - guess[k - 1]).norm();
  }
  if (!(cum.back() > 0.0)) {
    throw std::invalid_argument("optimizeInstanton: the two minima coincide");
  }
  const double width = betaHbar / (2.0 * options.betaHbarOmega);
  std::vector<VectorXd> x(static_cast<size_t>(P - 1));
  for (long j = 1; j < P; ++j) {
    const double tau = static_cast<double>(j) * dtau - 0.5 * betaHbar;
    const double f = 0.5 * (1.0 + std::tanh(tau / width));
    x[static_cast<size_t>(j - 1)] = alongPolyline(guess, cum, f);
  }

  // L-BFGS with a backtracking Armijo line search.
  ActionEval cur =
      evaluateAction(x, start, end, vEnds[0], vEnds[1], dtau, potential);
  std::deque<std::pair<std::vector<VectorXd>, std::vector<VectorXd>>> pairs;
  for (long it = 0; it < options.maxIterations; ++it) {
    inst.iterations = it;
    if (largestBeadNorm(cur.grad) / dtau < options.forceTolerance) {
      inst.converged = true;
      break;
    }
    std::vector<VectorXd> q = cur.grad;
    std::vector<double> alpha(pairs.size());
    for (size_t i = pairs.size(); i-- > 0;) {
      const double rho = 1.0 / dot(pairs[i].second, pairs[i].first);
      alpha[i] = rho * dot(pairs[i].first, q);
      for (size_t k = 0; k < q.size(); ++k) {
        q[k] -= alpha[i] * pairs[i].second[k];
      }
    }
    // Without history, half the inverse spring stiffness 2 / dtau.
    double gamma = dtau / 4.0;
    if (!pairs.empty()) {
      gamma = dot(pairs.back().first, pairs.back().second) /
              dot(pairs.back().second, pairs.back().second);
    }
    scale(q, gamma);
    for (size_t i = 0; i < pairs.size(); ++i) {
      const double rho = 1.0 / dot(pairs[i].second, pairs[i].first);
      const double beta = rho * dot(pairs[i].second, q);
      for (size_t k = 0; k < q.size(); ++k) {
        q[k] += (alpha[i] - beta) * pairs[i].first[k];
      }
    }
    // q is now the inverse-Hessian estimate times the gradient; step -q.
    double slope = -dot(cur.grad, q);
    if (!(slope < 0.0)) {
      pairs.clear();
      q = cur.grad;
      scale(q, dtau / 4.0);
      slope = -dot(cur.grad, q);
    }
    double step = 1.0;
    ActionEval next;
    std::vector<VectorXd> trial(x.size());
    bool accepted = false;
    for (int ls = 0; ls < 30; ++ls) {
      for (size_t k = 0; k < x.size(); ++k) {
        trial[k] = x[k] - step * q[k];
      }
      next = evaluateAction(trial, start, end, vEnds[0], vEnds[1], dtau,
                            potential);
      // Near the minimum the action changes by less than its round-off;
      // there a step that shrinks the gradient is progress too.
      const bool armijo = next.action <= cur.action + 1e-4 * step * slope;
      const bool flat = std::abs(next.action - cur.action) <=
                        1e-13 * std::max(1.0, std::abs(cur.action));
      if (armijo ||
          (flat && largestBeadNorm(next.grad) < largestBeadNorm(cur.grad))) {
        accepted = true;
        break;
      }
      step *= 0.5;
    }
    if (!accepted) {
      break;
    }
    std::vector<VectorXd> sk(x.size()), yk(x.size());
    for (size_t k = 0; k < x.size(); ++k) {
      sk[k] = trial[k] - x[k];
      yk[k] = next.grad[k] - cur.grad[k];
    }
    if (dot(sk, yk) > 0.0) {
      pairs.emplace_back(std::move(sk), std::move(yk));
      if (static_cast<long>(pairs.size()) > options.memory) {
        pairs.pop_front();
      }
    }
    x = std::move(trial);
    cur = std::move(next);
  }
  if (!inst.converged &&
      largestBeadNorm(cur.grad) / dtau < options.forceTolerance) {
    inst.converged = true;
  }

  inst.path.reserve(static_cast<size_t>(P + 1));
  inst.path.push_back(start);
  inst.path.insert(inst.path.end(), x.begin(), x.end());
  inst.path.push_back(end);
  inst.energies.reserve(static_cast<size_t>(P + 1));
  inst.energies.push_back(vEnds[0]);
  inst.energies.insert(inst.energies.end(), cur.v.begin(), cur.v.end());
  inst.energies.push_back(vEnds[1]);
  const double sWell = betaHbar * 0.5 * (vEnds[0] + vEnds[1]);
  inst.action = (cur.action - sWell) / kHbar;
  double s0 = 0.0;
  for (long j = 0; j < P; ++j) {
    s0 += (inst.path[static_cast<size_t>(j + 1)] -
           inst.path[static_cast<size_t>(j)])
              .squaredNorm();
  }
  inst.s0 = s0 / dtau;
  inst.symmetricEnough = std::abs(inst.asymmetry) * betaHbar / kHbar < 0.1;
  return inst;
}

void instantonSplitting(Instanton &inst, const BeadHessian &hessian,
                        const MatrixXd &hessStart, const MatrixXd &hessEnd) {
  const long P = static_cast<long>(inst.path.size()) - 1;
  if (P < 4 || !(inst.dtau > 0.0)) {
    throw std::invalid_argument("instantonSplitting: no optimised path");
  }
  const double dtau = inst.dtau;
  const double c = 1.0 / dtau;
  const long n = inst.path.front().size();
  const MatrixXd spring = 2.0 * c * MatrixXd::Identity(n, n);

  std::vector<MatrixXd> diag;
  diag.reserve(static_cast<size_t>(P - 1));
  for (long j = 1; j < P; ++j) {
    const MatrixXd h = hessian(j, inst.path[static_cast<size_t>(j)]);
    if (h.rows() != n || h.cols() != n) {
      throw std::runtime_error("instantonSplitting: bead Hessian size");
    }
    diag.push_back(spring + dtau * 0.5 * (h + h.transpose()));
  }
  const BlockChain chain(c, diag);

  auto wellLogDet = [&](const MatrixXd &h) {
    const std::vector<MatrixXd> d(static_cast<size_t>(P - 1),
                                  spring + dtau * 0.5 * (h + h.transpose()));
    const BlockChain well(c, d);
    if (well.sign() < 0) {
      throw std::runtime_error(
          "instantonSplitting: a well Hessian is not positive definite");
    }
    return well.logAbsDet();
  };
  const double logDetWell = 0.5 * (wellLogDet(hessStart) + wellLogDet(hessEnd));

  // The zero mode is the kink's translation in imaginary time, along the
  // discrete velocity v; det' J = det J (v^T J^-1 v) for v its eigenvector.
  std::vector<VectorXd> v(static_cast<size_t>(P - 1));
  for (long j = 1; j < P; ++j) {
    v[static_cast<size_t>(j - 1)] = inst.path[static_cast<size_t>(j + 1)] -
                                    inst.path[static_cast<size_t>(j - 1)];
  }
  scale(v, 1.0 / std::sqrt(dot(v, v)));
  const double vJv = dot(v, chain.solve(v));
  const int signPrime = chain.sign() * (vJv < 0.0 ? -1 : 1);
  if (signPrime < 0) {
    throw std::runtime_error(
        "instantonSplitting: the path is not a minimum of the action "
        "(a negative mode besides the kink's translation)");
  }
  inst.zeroMode = 1.0 / vJv;
  const double logDetPrime = chain.logAbsDet() + std::log(std::abs(vJv));

  // Next eigenvalue: inverse iteration orthogonal to v.
  std::vector<VectorXd> w(v.size());
  for (size_t k = 0; k < w.size(); ++k) {
    w[k].resize(n);
    for (long i = 0; i < n; ++i) {
      w[k](i) = std::sin(0.7 * static_cast<double>(k) +
                         1.3 * static_cast<double>(i) + 0.1);
    }
  }
  double lambda1 = 0.0;
  for (int it = 0; it < 40; ++it) {
    const double proj = dot(v, w);
    for (size_t k = 0; k < w.size(); ++k) {
      w[k] -= proj * v[k];
    }
    scale(w, 1.0 / std::sqrt(dot(w, w)));
    std::vector<VectorXd> z = chain.solve(w);
    lambda1 = 1.0 / dot(w, z);
    w = std::move(z);
  }
  inst.modeSeparation = std::abs(lambda1 / inst.zeroMode);

  inst.delta0 = 2.0 * kHbar *
                std::sqrt(inst.s0 / (2.0 * std::numbers::pi * kHbar * dtau)) *
                std::exp(0.5 * (logDetWell - logDetPrime) - inst.action);
}

double crossoverTemperature(const MatrixXd &hessSaddle) {
  const auto es = symEigen(hessSaddle);
  const double lambda = es.eigenvalues()(0);
  if (!(lambda < 0.0)) {
    throw std::invalid_argument(
        "crossoverTemperature: the saddle Hessian has no negative eigenvalue");
  }
  return kHbar * std::sqrt(-lambda) / (2.0 * std::numbers::pi * kBoltzmann);
}

namespace {

using Ring = std::vector<VectorXd>;

Ring zeroRing(long n, long f) {
  return Ring(static_cast<size_t>(n), VectorXd::Zero(f));
}

/// J v for the ring Hessian of bead blocks `h` and spring constant c.
Ring ringApply(const std::vector<MatrixXd> &h, double c, const Ring &v) {
  const size_t n = v.size();
  Ring out(n);
  for (size_t j = 0; j < n; ++j) {
    out[j] =
        h[j] * v[j] + c * (2.0 * v[j] - v[(j + n - 1) % n] - v[(j + 1) % n]);
  }
  return out;
}

/// Block LU of the open chain T: blocks h_j + 2 c I on the diagonal, -c I
/// between neighbours, no closure. Solves and, on request, the inertia and
/// log-determinant from the Schur complements (Haynsworth: the inertia of T
/// is the sum over its Schur blocks).
class OpenChain {
public:
  OpenChain(double c, const std::vector<MatrixXd> &h, bool spectrum)
      : c_(c) {
    const size_t n = h.size();
    lu_.reserve(n);
    for (size_t k = 0; k < n; ++k) {
      const long f = h[k].rows();
      MatrixXd d =
          0.5 * (h[k] + h[k].transpose()) + 2.0 * c * MatrixXd::Identity(f, f);
      if (k > 0) {
        // S_k = D_k - c^2 S_{k-1}^{-1}, symmetric like its predecessor.
        const MatrixXd inv = lu_.back().inverse();
        d -= c * c * 0.5 * (inv + inv.transpose());
      }
      if (spectrum) {
        const auto es = symEigen(d, false);
        for (long i = 0; i < f; ++i) {
          const double lam = es.eigenvalues()(i);
          if (lam == 0.0) {
            throw std::runtime_error("instanton: singular chain Hessian block");
          }
          if (lam < 0.0) {
            ++negative_;
            sign_ = -sign_;
          }
          logAbsDet_ += std::log(std::abs(lam));
        }
      }
      lu_.emplace_back(d);
    }
  }
  double logAbsDet() const { return logAbsDet_; }
  int sign() const { return sign_; }
  long negative() const { return negative_; }
  Ring solve(const Ring &b) const {
    const size_t m = b.size();
    Ring y(m), x(m);
    y[0] = b[0];
    for (size_t k = 1; k < m; ++k) {
      y[k] = b[k] + c_ * lu_[k - 1].solve(y[k - 1]);
    }
    x[m - 1] = lu_[m - 1].solve(y[m - 1]);
    for (size_t k = m - 1; k-- > 0;) {
      x[k] = lu_[k].solve(y[k] + c_ * x[k + 1]);
    }
    return x;
  }

private:
  double c_;
  std::vector<Eigen::PartialPivLU<ColMajorXd>> lu_;
  double logAbsDet_ = 0.0;
  int sign_ = 1;
  long negative_ = 0;
};

/// The ring Hessian J = T + G K G^T: the closure blocks (-c I between bead
/// 0 and bead N - 1) and extra symmetric rank-one terms kappa_i u_i u_i^T,
/// through the open chain and a small dense matrix of size 2 f + extras.
/// Woodbury gives J^{-1} b, the determinant lemma ln |det J| and Haynsworth
/// the inertia, all O(N f^3).
class ClosedRing {
public:
  ClosedRing(double c, const std::vector<MatrixXd> &h,
             const std::vector<Ring> &extras, const std::vector<double> &kappas,
             bool spectrum)
      : c_(c),
        n_(static_cast<long>(h.size())),
        f_(h.front().rows()),
        chain_(c, h, spectrum),
        extras_(extras),
        kappas_(kappas) {
    const long m = 2 * f_ + static_cast<long>(extras.size());
    // K^{-1}: the closure block inverts to [[0, -I/c], [-I/c, 0]].
    kinv_ = MatrixXd::Zero(m, m);
    kinv_.block(0, f_, f_, f_) = -MatrixXd::Identity(f_, f_) / c;
    kinv_.block(f_, 0, f_, f_) = -MatrixXd::Identity(f_, f_) / c;
    for (size_t i = 0; i < extras.size(); ++i) {
      kinv_(2 * f_ + static_cast<long>(i), 2 * f_ + static_cast<long>(i)) =
          1.0 / kappas[i];
    }
    // G^T T^{-1} G, column by column: the pieces of T^{-1} g at bead 0,
    // bead N - 1 and along every extra vector.
    gtg_ = MatrixXd::Zero(m, m);
    for (long col = 0; col < m; ++col) {
      const Ring x = chain_.solve(column(col));
      gtg_.col(col) = pieces(x);
    }
    // (K^{-1} + G^T T^{-1} G) for Woodbury; I + K G^T T^{-1} G for the
    // determinant.
    woodbury_.compute(ColMajorXd(kinv_ + gtg_));
    if (spectrum) {
      MatrixXd k = MatrixXd::Zero(m, m);
      k.block(0, f_, f_, f_) = -c * MatrixXd::Identity(f_, f_);
      k.block(f_, 0, f_, f_) = -c * MatrixXd::Identity(f_, f_);
      for (size_t i = 0; i < extras.size(); ++i) {
        k(2 * f_ + static_cast<long>(i), 2 * f_ + static_cast<long>(i)) =
            kappas[i];
      }
      const Eigen::PartialPivLU<ColMajorXd> lu(MatrixXd::Identity(m, m) +
                                               k * gtg_);
      double logDet = 0.0;
      int sign = static_cast<int>(std::lround(lu.permutationP().determinant()));
      for (long i = 0; i < m; ++i) {
        const double u = lu.matrixLU()(i, i);
        if (u == 0.0) {
          throw std::runtime_error("instanton: singular ring closure");
        }
        sign *= u < 0.0 ? -1 : 1;
        logDet += std::log(std::abs(u));
      }
      logAbsDet_ = chain_.logAbsDet() + logDet;
      sign_ = chain_.sign() * sign;
      // neg(J) = neg(T) + neg(S) - neg(-K^{-1}), S = -K^{-1} - G^T T^{-1} G.
      const MatrixXd s = -kinv_ - gtg_;
      const auto es = symEigen(s, false);
      const auto ek = symEigen(-kinv_, false);
      negative_ = chain_.negative() + (es.eigenvalues().array() < 0.0).count() -
                  (ek.eigenvalues().array() < 0.0).count();
    }
  }
  Ring solve(const Ring &b) const {
    const Ring tb = chain_.solve(b);
    const VectorXd y = woodbury_.solve(pieces(tb));
    const Ring gy = chain_.solve(expand(y));
    Ring x(tb.size());
    for (size_t j = 0; j < x.size(); ++j) {
      x[j] = tb[j] - gy[j];
    }
    return x;
  }
  double logAbsDet() const { return logAbsDet_; }
  int sign() const { return sign_; }
  long negative() const { return negative_; }

private:
  /// Column `col` of G as a ring vector.
  Ring column(long col) const {
    if (col < 2 * f_) {
      Ring g = zeroRing(n_, f_);
      g[static_cast<size_t>(col < f_ ? 0 : n_ - 1)](col % f_) = 1.0;
      return g;
    }
    return extras_[static_cast<size_t>(col - 2 * f_)];
  }
  /// G^T x.
  VectorXd pieces(const Ring &x) const {
    VectorXd out(2 * f_ + static_cast<long>(extras_.size()));
    out.head(f_) = x[0];
    out.segment(f_, f_) = x[static_cast<size_t>(n_ - 1)];
    for (size_t i = 0; i < extras_.size(); ++i) {
      out(2 * f_ + static_cast<long>(i)) = dot(extras_[i], x);
    }
    return out;
  }
  /// G y.
  Ring expand(const VectorXd &y) const {
    Ring g = zeroRing(n_, f_);
    g[0] += y.head(f_);
    g[static_cast<size_t>(n_ - 1)] += y.segment(f_, f_);
    for (size_t i = 0; i < extras_.size(); ++i) {
      const double w = y(2 * f_ + static_cast<long>(i));
      for (long j = 0; j < n_; ++j) {
        g[static_cast<size_t>(j)] += w * extras_[i][static_cast<size_t>(j)];
      }
    }
    return g;
  }

  double c_;
  long n_, f_;
  OpenChain chain_;
  std::vector<Ring> extras_;
  std::vector<double> kappas_;
  MatrixXd kinv_, gtg_;
  Eigen::PartialPivLU<ColMajorXd> woodbury_;
  double logAbsDet_ = 0.0;
  int sign_ = 1;
  long negative_ = 0;
};

/// The normalised imaginary-time translation, tau_j = (q_{j+1} - q_{j-1}) / 2.
Ring timeTranslation(const Ring &x) {
  const size_t n = x.size();
  Ring tau(n);
  for (size_t j = 0; j < n; ++j) {
    tau[j] = 0.5 * (x[(j + 1) % n] - x[(j + n - 1) % n]);
  }
  const double norm = std::sqrt(dot(tau, tau));
  if (!(norm > 0.0)) {
    throw std::runtime_error("instanton: the ring has collapsed to a point");
  }
  scale(tau, 1.0 / norm);
  return tau;
}

struct RitzPair {
  double theta = 0.0;
  Ring vector;
};

/// The lowest `count` Ritz pairs of the ring Hessian by Lanczos with full
/// reorthogonalisation, started from `start`, ascending.
std::vector<RitzPair> lowestRingModes(const std::vector<MatrixXd> &h, double c,
                                      const Ring &start, long steps,
                                      long count) {
  const size_t n = start.size();
  std::vector<Ring> basis;
  std::vector<double> alpha, beta;
  Ring q = start;
  scale(q, 1.0 / std::sqrt(dot(q, q)));
  for (long k = 0; k < steps; ++k) {
    basis.push_back(q);
    Ring w = ringApply(h, c, q);
    alpha.push_back(dot(w, q));
    for (int pass = 0; pass < 2; ++pass) {
      for (const auto &b : basis) {
        const double p = dot(w, b);
        for (size_t j = 0; j < n; ++j) {
          w[j] -= p * b[j];
        }
      }
    }
    const double bnorm = std::sqrt(dot(w, w));
    if (!(bnorm > 1e-10) || k + 1 == steps) {
      break;
    }
    beta.push_back(bnorm);
    scale(w, 1.0 / bnorm);
    q = std::move(w);
  }
  const long m = static_cast<long>(alpha.size());
  MatrixXd t = MatrixXd::Zero(m, m);
  for (long i = 0; i < m; ++i) {
    t(i, i) = alpha[static_cast<size_t>(i)];
    if (i + 1 < m) {
      t(i, i + 1) = t(i + 1, i) = beta[static_cast<size_t>(i)];
    }
  }
  const auto es = symEigen(t);
  std::vector<RitzPair> out;
  for (long r = 0; r < std::min(count, m); ++r) {
    const VectorXd y = es.eigenvectors().col(r);
    RitzPair pair;
    pair.theta = es.eigenvalues()(r);
    pair.vector = zeroRing(static_cast<long>(n), start.front().size());
    for (long i = 0; i < m; ++i) {
      for (size_t j = 0; j < n; ++j) {
        pair.vector[j] += y(i) * basis[static_cast<size_t>(i)][j];
      }
    }
    scale(pair.vector, 1.0 / std::sqrt(dot(pair.vector, pair.vector)));
    out.push_back(std::move(pair));
  }
  return out;
}

/// The lowest eigenpair alone; `mode` returns holding the eigenvector.
double lowestRingMode(const std::vector<MatrixXd> &h, double c, Ring &mode,
                      long steps) {
  auto pairs = lowestRingModes(h, c, mode, steps, 1);
  mode = std::move(pairs.front().vector);
  return pairs.front().theta;
}

/// Bofill's update of a bead Hessian from a step and its gradient change.
void bofillUpdate(MatrixXd &h, const VectorXd &dq, const VectorXd &dg) {
  const double dq2 = dq.squaredNorm();
  if (!(dq2 > 1e-24)) {
    return;
  }
  const VectorXd r = dg - h * dq;
  const double r2 = r.squaredNorm();
  const double rq = r.dot(dq);
  if (!(r2 > 0.0)) {
    return;
  }
  const double phi = rq * rq / (r2 * dq2);
  MatrixXd update = (r * dq.transpose() + dq * r.transpose()) / dq2 -
                    rq * dq * dq.transpose() / (dq2 * dq2);
  if (std::abs(rq) > 1e-12 * std::sqrt(r2 * dq2)) {
    update = phi * (r * r.transpose()) / rq + (1.0 - phi) * update;
  }
  h += update;
  h = 0.5 * (h + h.transpose()).eval();
}

struct RingEval {
  double u = 0.0; // U_N
  std::vector<double> v;
  Ring gradV; // dV/dq at each bead
  Ring grad;  // dU_N/dq at each bead
};

/// U_N, V and both gradients over the ring. With `symmetric` the ring is
/// taken as invariant under imaginary-time reversal, bead j the mirror of
/// bead N - j, so only beads 0..N/2 go to the potential.
RingEval evaluateRing(const Ring &x, double c, const BatchPotential &potential,
                      bool symmetric) {
  RingEval out;
  const size_t n = x.size();
  if (symmetric && n % 2 == 0) {
    const size_t half = n / 2;
    const Ring xs(x.begin(), x.begin() + static_cast<long>(half) + 1);
    std::vector<double> vh;
    Ring gh;
    potential(xs, vh, gh);
    if (vh.size() != half + 1 || gh.size() != half + 1) {
      throw std::runtime_error(
          "rate instanton: potential returned the wrong count");
    }
    out.v.resize(n);
    out.gradV.resize(n);
    for (size_t j = 0; j < n; ++j) {
      const size_t k = j <= half ? j : n - j;
      out.v[j] = vh[k];
      out.gradV[j] = gh[k];
    }
  } else {
    potential(x, out.v, out.gradV);
    if (out.v.size() != n || out.gradV.size() != n) {
      throw std::runtime_error(
          "rate instanton: potential returned the wrong count");
    }
  }
  out.grad.resize(n);
  double spring = 0.0;
  for (size_t j = 0; j < n; ++j) {
    const VectorXd &prev = x[(j + n - 1) % n];
    const VectorXd &next = x[(j + 1) % n];
    out.grad[j] = out.gradV[j] + c * (2.0 * x[j] - prev - next);
    spring += (next - x[j]).squaredNorm();
    out.u += out.v[j];
  }
  out.u += 0.5 * c * spring;
  return out;
}

/// Averages bead j with its mirror N - j.
void symmetrise(Ring &x) {
  const size_t n = x.size();
  for (size_t j = 1; j < n - j; ++j) {
    x[j] = 0.5 * (x[j] + x[n - j]);
    x[n - j] = x[j];
  }
}

double ringLength2(const Ring &x) {
  double b = 0.0;
  for (size_t j = 0; j < x.size(); ++j) {
    b += (x[(j + 1) % x.size()] - x[j]).squaredNorm();
  }
  return b;
}

/// Distance along sign * dir from the saddle at which V has dropped by
/// `drop`, from a scan in steps of h; the lowest point found when V never
/// drops that far before rising again.
double turningDistance(const VectorXd &saddle, const VectorXd &dir,
                       double vSaddle, double drop, double sign, double h,
                       long maxPoints, const BatchPotential &potential) {
  Ring pts;
  for (long k = 1; k <= maxPoints; ++k) {
    pts.push_back(saddle + sign * h * static_cast<double>(k) * dir);
  }
  std::vector<double> v;
  Ring g;
  potential(pts, v, g);
  double prevV = vSaddle, prevS = 0.0, best = 0.0, bestV = vSaddle;
  for (long k = 0; k < maxPoints; ++k) {
    const double sk = h * static_cast<double>(k + 1);
    const double vk = v[static_cast<size_t>(k)];
    if (vk <= vSaddle - drop) {
      const double t = (prevV - (vSaddle - drop)) / (prevV - vk);
      return prevS + t * (sk - prevS);
    }
    if (vk < bestV) {
      bestV = vk;
      best = sk;
    }
    if (vk > prevV + 1e-12 && k > 0) {
      break;
    }
    prevV = vk;
    prevS = sk;
  }
  return best;
}

double sideDrop(const VectorXd &saddle, const VectorXd &dir, double vSaddle,
                double sign, double h, long maxPoints,
                const BatchPotential &potential) {
  Ring pts;
  for (long k = 1; k <= maxPoints; ++k) {
    pts.push_back(saddle + sign * h * static_cast<double>(k) * dir);
  }
  std::vector<double> v;
  Ring g;
  potential(pts, v, g);
  double lowest = vSaddle;
  for (long k = 0; k < maxPoints; ++k) {
    const double vk = v[static_cast<size_t>(k)];
    if (vk > lowest + 1e-12 && k > 0 && lowest < vSaddle) {
      return vSaddle - lowest;
    }
    lowest = std::min(lowest, vk);
  }
  return vSaddle - lowest; // still falling: the drop to the last point
}

/// Cumulative arc length of a polyline.
std::vector<double> arcLengths(const Ring &path) {
  std::vector<double> s(path.size(), 0.0);
  for (size_t k = 1; k < path.size(); ++k) {
    s[k] = s[k - 1] + (path[k] - path[k - 1]).norm();
  }
  return s;
}

/// The point of a polyline at arc length `target`.
VectorXd atArcLength(const Ring &path, const std::vector<double> &s,
                     double target) {
  const auto it = std::upper_bound(s.begin(), s.end(), target);
  const size_t k = std::clamp<size_t>(
      static_cast<size_t>(std::distance(s.begin(), it)), 1, path.size() - 1);
  const double seg = s[k] - s[k - 1];
  const double t = seg > 0.0 ? (target - s[k - 1]) / seg : 0.0;
  return path[k - 1] + std::clamp(t, 0.0, 1.0) * (path[k] - path[k - 1]);
}

/// Turning points of the classical orbit at energy E in the inverted
/// potential: the crossings V(s) = E nearest the barrier top on each side.
std::pair<double, double> turningPoints(const Profile &p, double sTop,
                                        double e) {
  const double s0 = p.s().front(), s1 = p.s().back();
  auto cross = [&](double from, double to) {
    // V(from) > e >= V(to) somewhere between: scan, then bisect.
    const int steps = 2000;
    double a = from, b = to;
    for (int k = 1; k <= steps; ++k) {
      const double sk = from + (to - from) * static_cast<double>(k) / steps;
      if (p(sk) <= e) {
        a = from + (to - from) * static_cast<double>(k - 1) / steps;
        b = sk;
        break;
      }
      if (k == steps) {
        return to;
      }
    }
    for (int k = 0; k < 60; ++k) {
      const double m = 0.5 * (a + b);
      (p(m) > e ? a : b) = m;
    }
    return 0.5 * (a + b);
  };
  return {cross(sTop, s0), cross(sTop, s1)};
}

/// Half period int_{s-}^{s+} ds / sqrt(2 (V(s) - E)) with the inverse square
/// roots at the turning points removed by s = mid - half cos(phi).
double halfPeriod(const Profile &p, double sMinus, double sPlus, double e,
                  std::vector<double> *cumulative = nullptr,
                  std::vector<double> *positions = nullptr) {
  const int m = 4000;
  const double mid = 0.5 * (sMinus + sPlus), half = 0.5 * (sPlus - sMinus);
  double total = 0.0;
  if (cumulative) {
    cumulative->assign(1, 0.0);
    positions->assign(1, sMinus);
  }
  for (int k = 0; k < m; ++k) {
    const double phi = std::numbers::pi * (static_cast<double>(k) + 0.5) / m;
    const double s = mid - half * std::cos(phi);
    const double under = 2.0 * (p(s) - e);
    const double integrand =
        half * std::sin(phi) / std::sqrt(std::max(under, 1e-300));
    total += integrand * std::numbers::pi / m;
    if (cumulative) {
      cumulative->push_back(total);
      positions->push_back(mid -
                           half * std::cos(std::numbers::pi * (k + 1.0) / m));
    }
  }
  return total;
}

} // namespace

RingSpectrum ringSpectrum(const std::vector<MatrixXd> &beadHessians, double c,
                          const std::vector<VectorXd> &tau) {
  if (beadHessians.empty() || beadHessians.size() != tau.size()) {
    throw std::invalid_argument(
        "ringSpectrum: N bead Hessians and N tau blocks");
  }
  const ClosedRing ring(c, beadHessians, {tau}, {1.0}, true);
  RingSpectrum out;
  out.logDetPrime = ring.logAbsDet();
  out.signDetPrime = ring.sign();
  out.negativeModes = ring.negative();
  out.zeroEigenvalue = dot(tau, ringApply(beadHessians, c, tau));
  return out;
}

std::vector<VectorXd> ringFromPath(const std::vector<VectorXd> &path,
                                   const std::vector<double> &energies,
                                   double betaHbar, long beads) {
  if (path.size() < 3 || path.size() != energies.size() || beads < 4 ||
      !(betaHbar > 0.0)) {
    throw std::invalid_argument("ringFromPath: a path of at least three points "
                                "with energies, N >= 4 and "
                                "beta hbar > 0");
  }
  const std::vector<double> s = arcLengths(path);
  const Profile p(s, energies);
  // The barrier top on a fine grid.
  double sTop = s.front(), vTop = -std::numeric_limits<double>::infinity();
  const int grid = 4000;
  for (int k = 0; k <= grid; ++k) {
    const double sk =
        s.front() + (s.back() - s.front()) * static_cast<double>(k) / grid;
    if (p(sk) > vTop) {
      vTop = p(sk);
      sTop = sk;
    }
  }
  const double vLow = std::max(energies.front(), energies.back());
  if (!(vTop > vLow)) {
    throw std::invalid_argument("ringFromPath: the path has no barrier");
  }
  auto period = [&](double e) {
    const auto [sm, sp] = turningPoints(p, sTop, e);
    return 2.0 * halfPeriod(p, sm, sp, e);
  };
  // The crossover along the path from a parabola through the three points
  // around the barrier top of the input, not of the interpolant.
  {
    size_t top = 0;
    for (size_t k = 1; k < energies.size(); ++k) {
      if (energies[k] > energies[top]) {
        top = k;
      }
    }
    if (top == 0 || top + 1 == energies.size()) {
      throw std::invalid_argument(
          "ringFromPath: the barrier top is an end of the path");
    }
    const double h1 = s[top] - s[top - 1], h2 = s[top + 1] - s[top];
    const double curvature =
        2.0 *
        (h1 * energies[top + 1] - (h1 + h2) * energies[top] +
         h2 * energies[top - 1]) /
        (h1 * h2 * (h1 + h2));
    if (!(curvature < 0.0)) {
      throw std::invalid_argument("ringFromPath: no curvature at the top");
    }
    const double tc = kHbar * std::sqrt(-curvature) / (2.0 * std::numbers::pi);
    if (!(kHbar / betaHbar < tc)) {
      throw std::invalid_argument(
          "ringFromPath: the temperature is at or above the crossover along "
          "this path");
    }
  }
  // Bracket the orbit energy geometrically above the lower end: the period
  // grows only logarithmically as E approaches a well bottom.
  double eHi = vTop - 1e-9 * (vTop - vLow);
  double eLo = vLow + 1e-14 * (vTop - vLow);
  if (period(eLo) < betaHbar) {
    // The path does not reach far enough down for this temperature; the
    // lowest orbit it holds is the best start.
    eHi = eLo;
  }
  for (int k = 0; k < 200 && eHi > eLo; ++k) {
    const double e = vLow + std::sqrt((eLo - vLow) * (eHi - vLow));
    (period(e) > betaHbar ? eLo : eHi) = e;
    if (eHi - eLo < 1e-15 * (vTop - vLow)) {
      break;
    }
  }
  const double e = 0.5 * (eLo + eHi);
  const auto [sm, sp] = turningPoints(p, sTop, e);
  std::vector<double> tau, pos;
  const double half = halfPeriod(p, sm, sp, e, &tau, &pos);
  // Bead j sits at imaginary time j beta hbar / N along the half orbit from
  // the reactant side turning point (j = 0) to the product side (j = N / 2),
  // and bead N - j mirrors bead j.
  std::vector<VectorXd> ring(static_cast<size_t>(beads));
  for (long j = 0; j <= beads / 2; ++j) {
    const double t =
        std::min(half, half * 2.0 * static_cast<double>(j) / beads);
    const auto it = std::upper_bound(tau.begin(), tau.end(), t);
    const size_t k = std::clamp<size_t>(
        static_cast<size_t>(std::distance(tau.begin(), it)), 1, tau.size() - 1);
    const double seg = tau[k] - tau[k - 1];
    const double w = seg > 0.0 ? (t - tau[k - 1]) / seg : 0.0;
    const double sj =
        pos[k - 1] + std::clamp(w, 0.0, 1.0) * (pos[k] - pos[k - 1]);
    ring[static_cast<size_t>(j)] = atArcLength(path, s, sj);
    if (j > 0 && j < beads - j) {
      ring[static_cast<size_t>(beads - j)] = ring[static_cast<size_t>(j)];
    }
  }
  return ring;
}

double wkbLogRateAlongPath(const Profile &profile, double beta,
                           double hwReactant) {
  if (!(beta > 0.0) || !(hwReactant > 0.0)) {
    throw std::invalid_argument("wkbLogRateAlongPath: beta and hbar omega > 0");
  }
  const double vR = profile.v().front();
  double vTop = vR;
  const int grid = 2000;
  for (int k = 0; k <= grid; ++k) {
    const double sk =
        profile.s().front() + (profile.s().back() - profile.s().front()) *
                                  static_cast<double>(k) / grid;
    vTop = std::max(vTop, profile(sk));
  }
  const double barrier = vTop - vR;
  if (!(barrier > 0.0)) {
    throw std::invalid_argument(
        "wkbLogRateAlongPath: no barrier above the reactant");
  }
  // Curvature at the top from the profile, for Kemble's continuation above
  // the barrier: theta = -pi (E - V_top) / (hbar omega_b).
  double sTop = profile.s().front();
  for (int k = 0; k <= grid; ++k) {
    const double sk =
        profile.s().front() + (profile.s().back() - profile.s().front()) *
                                  static_cast<double>(k) / grid;
    if (profile(sk) >= vTop) {
      sTop = sk;
    }
  }
  const double ds = 1e-3 * (profile.s().back() - profile.s().front());
  const double curvature =
      std::max(1e-12, -(profile(sTop + ds) - 2.0 * vTop + profile(sTop - ds)) /
                          (ds * ds));
  const double hwB = kHbar * std::sqrt(curvature);
  // int P(E) exp(-beta (E - V_R)) dE from the reactant minimum up to where
  // the Boltzmann factor has died, on a grid dense below the barrier.
  const double eMax = barrier + 40.0 / beta;
  const int points = 600;
  double logTerms = -std::numeric_limits<double>::infinity();
  auto logAdd = [](double a, double b) {
    if (a == -std::numeric_limits<double>::infinity()) {
      return b;
    }
    const double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
  };
  double prevLog = -std::numeric_limits<double>::infinity(), prevE = 0.0;
  for (int k = 0; k <= points; ++k) {
    const double e = eMax * static_cast<double>(k) / points;
    double theta;
    if (e < barrier) {
      theta = wkbAction(profile, vR + e);
    } else {
      theta = -std::numbers::pi * (e - barrier) / hwB;
    }
    // ln P = -ln(1 + exp(2 theta)).
    const double logP =
        theta > 20.0 ? -2.0 * theta : -std::log1p(std::exp(2.0 * theta));
    const double logF = logP - beta * e;
    if (k > 0) {
      // trapezoid in log space
      const double segment =
          std::log(0.5 * (e - prevE)) + logAdd(prevLog, logF);
      logTerms = logAdd(logTerms, segment);
    }
    prevLog = logF;
    prevE = e;
  }
  const double logFlux = logTerms - std::log(2.0 * std::numbers::pi * kHbar);
  return logFlux + std::log(2.0 * std::sinh(0.5 * beta * hwReactant));
}

RateInstanton optimizeRateInstanton(const VectorXd &saddle,
                                    const MatrixXd &hessSaddle, double beta,
                                    std::vector<VectorXd> guess,
                                    const BatchPotential &potential,
                                    const RateInstantonOptions &options,
                                    std::vector<MatrixXd> beadHessians) {
  const long N = options.beads;
  if (N < 4 || !(beta > 0.0) || hessSaddle.rows() != saddle.size()) {
    throw std::invalid_argument(
        "optimizeRateInstanton: need N >= 4, beta > 0 and a saddle Hessian of "
        "the saddle's dimension");
  }
  RateInstanton inst;
  inst.beta = beta;
  inst.betaN = beta / static_cast<double>(N);
  inst.temperature = 1.0 / (kBoltzmann * beta);
  inst.crossover = crossoverTemperature(hessSaddle);
  if (!(inst.temperature < inst.crossover)) {
    throw std::invalid_argument(
        "optimizeRateInstanton: T is at or above the crossover temperature; "
        "the ring collapses onto the saddle and classical transition-state "
        "theory applies");
  }
  const double bnh = inst.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);
  const long f = saddle.size();
  const MatrixXd hs = 0.5 * (hessSaddle + hessSaddle.transpose());
  const auto es = symEigen(hs);
  const VectorXd dir = es.eigenvectors().col(0);

  if (static_cast<long>(guess.size()) != N) {
    std::vector<double> v0;
    Ring g0;
    potential({saddle}, v0, g0);
    const double vS = v0.at(0);
    // Scan in steps of a quarter of the length over which the barrier's
    // curvature drops V by kB T_c.
    const double h = 0.25 * std::sqrt(2.0 * kBoltzmann * inst.crossover /
                                      -es.eigenvalues()(0));
    const long pts = 200;
    const double dPlus = sideDrop(saddle, dir, vS, 1.0, h, pts, potential);
    const double dMinus = sideDrop(saddle, dir, vS, -1.0, h, pts, potential);
    const double dMin = std::min(dPlus, dMinus);
    const double drop = (1.0 - inst.temperature / inst.crossover) *
                        (dMin > 0.0 ? dMin : kBoltzmann * inst.crossover);
    const double sPlus =
        turningDistance(saddle, dir, vS, drop, 1.0, h, pts, potential);
    const double sMinus =
        turningDistance(saddle, dir, vS, drop, -1.0, h, pts, potential);
    guess.resize(static_cast<size_t>(N));
    for (long j = 0; j < N; ++j) {
      const double ct =
          std::cos(2.0 * std::numbers::pi * static_cast<double>(j) /
                   static_cast<double>(N));
      guess[static_cast<size_t>(j)] =
          saddle + dir * (ct >= 0.0 ? sPlus * ct : sMinus * ct);
    }
  }
  if (static_cast<long>(beadHessians.size()) != N) {
    beadHessians.assign(static_cast<size_t>(N), hs);
  }
  for (auto &h : beadHessians) {
    if (h.rows() != f || h.cols() != f) {
      throw std::invalid_argument(
          "optimizeRateInstanton: a bead Hessian of the wrong size");
    }
  }

  Ring x = std::move(guess);
  const bool symmetric = options.timeReversalSymmetric && N % 2 == 0;
  if (symmetric) {
    symmetrise(x);
    for (long j = 1; j < N - j; ++j) {
      beadHessians[static_cast<size_t>(N - j)] =
          beadHessians[static_cast<size_t>(j)];
    }
  }
  RingEval cur = evaluateRing(x, c, potential, symmetric);
  Ring mode(x.size(), dir);
  double trust = options.trustRadius;
  double lowest = 0.0;
  // A first-order saddle has one negative mode over the zero mode. Ritz
  // pairs that overlap the imaginary-time translation are that zero mode.
  auto isTranslation = [](const Ring &v, const Ring &tau) {
    return std::abs(dot(v, tau)) > 0.9;
  };
  const long ritzCount = 6;
  for (long it = 0; it < options.maxIterations; ++it) {
    inst.iterations = it;
    const Ring tau = timeTranslation(x);
    std::vector<RitzPair> ritz =
        lowestRingModes(beadHessians, c, mode, options.lanczosSteps, ritzCount);
    mode = ritz.front().vector;
    lowest = ritz.front().theta;
    bool extraNegative = false;
    for (size_t r = 1; r < ritz.size(); ++r) {
      if (ritz[r].theta < -1e-8 && !isTranslation(ritz[r].vector, tau)) {
        extraNegative = true;
      }
    }
    if (lowest < 0.0 && !extraNegative &&
        largestBeadNorm(cur.grad) < options.forceTolerance) {
      inst.converged = true;
      break;
    }
    // Eigenvector following: the Newton step climbs along the lowest mode
    // on its own once that mode is negative; while it is still positive its
    // sign is flipped, and every further negative mode (the zero mode aside)
    // is flipped too, so the step descends along it towards a first-order
    // saddle. The imaginary-time translation is held with a spring-sized
    // curvature and its component removed from the step.
    std::vector<Ring> extras{tau};
    std::vector<double> kappas{c};
    if (lowest > 0.0) {
      extras.push_back(mode);
      kappas.push_back(-2.0 * lowest);
    }
    for (size_t r = 1; r < ritz.size(); ++r) {
      if (ritz[r].theta < -1e-8 && !isTranslation(ritz[r].vector, tau)) {
        extras.push_back(ritz[r].vector);
        kappas.push_back(-2.0 * ritz[r].theta);
      }
    }
    const ClosedRing ring(c, beadHessians, extras, kappas, false);
    Ring step = ring.solve(cur.grad);
    scale(step, -1.0);
    const double along = dot(step, tau);
    for (size_t j = 0; j < step.size(); ++j) {
      step[j] -= along * tau[j];
    }
    const double big = largestBeadNorm(step);
    if (big > trust) {
      scale(step, trust / big);
    }
    const Ring js = ringApply(beadHessians, c, step);
    const double predicted = dot(cur.grad, step) + 0.5 * dot(step, js);
    Ring trial(x.size());
    for (size_t j = 0; j < x.size(); ++j) {
      trial[j] = x[j] + step[j];
    }
    if (symmetric) {
      symmetrise(trial);
      for (size_t j = 0; j < x.size(); ++j) {
        step[j] = trial[j] - x[j];
      }
    }
    RingEval next = evaluateRing(trial, c, potential, symmetric);
    if (options.updateHessians) {
      for (size_t j = 0; j < x.size(); ++j) {
        bofillUpdate(beadHessians[j], step[j], next.gradV[j] - cur.gradV[j]);
      }
    }
    const double ratio =
        std::abs(predicted) > 1e-30 ? (next.u - cur.u) / predicted : 1.0;
    if (ratio > 0.75 && ratio < 1.25 && big >= 0.99 * trust) {
      trust = std::min(2.0 * trust, options.maxTrustRadius);
    } else if (ratio < 0.25 || ratio > 1.75) {
      trust = std::max(0.5 * trust, 1e-4);
    }
    x = std::move(trial);
    cur = std::move(next);
  }
  if (!inst.converged) {
    lowest = lowestRingMode(beadHessians, c, mode, options.lanczosSteps);
    inst.converged =
        lowest < 0.0 && largestBeadNorm(cur.grad) < options.forceTolerance;
    if (inst.converged) {
      inst.iterations = options.maxIterations;
    }
  }
  inst.beads = x;
  inst.energies = cur.v;
  inst.gradients = cur.gradV;
  inst.hessians = std::move(beadHessians);
  inst.ringPotential = cur.u;
  inst.bN = ringLength2(x);
  inst.lowestEigenvalue = lowest;
  return inst;
}

void instantonRate(RateInstanton &inst, const RingBeadHessian &hessian,
                   const MatrixXd &hessReactant, double vReactant,
                   const MatrixXd &hessSaddle, double vSaddle,
                   const MatrixXd &rigidBasis) {
  const long N = static_cast<long>(inst.beads.size());
  if (N < 4 || !(inst.betaN > 0.0)) {
    throw std::invalid_argument("instantonRate: no optimised ring");
  }
  const long f = inst.beads.front().size();
  const long m = rigidBasis.size() > 0 ? rigidBasis.cols() : 0;
  if (m > 0 && rigidBasis.rows() != f) {
    throw std::invalid_argument(
        "instantonRate: rigid basis of the wrong dimension");
  }
  const double bnh = inst.betaN * kHbar;
  const double c = 1.0 / (bnh * bnh);
  std::vector<MatrixXd> h(static_cast<size_t>(N));
  for (long j = 0; j < N; ++j) {
    MatrixXd hj = hessian(j, inst.beads[static_cast<size_t>(j)]);
    if (hj.rows() != f || hj.cols() != f) {
      throw std::runtime_error("instantonRate: bead Hessian size");
    }
    h[static_cast<size_t>(j)] = 0.5 * (hj + hj.transpose());
  }
  // det' leaves out the imaginary-time translation and the rigid modes,
  // each a null vector of J: uniform rigid displacements of every bead.
  std::vector<Ring> dropped{timeTranslation(inst.beads)};
  std::vector<double> kappas{1.0};
  for (long i = 0; i < m; ++i) {
    dropped.emplace_back(
        static_cast<size_t>(N),
        (rigidBasis.col(i) / std::sqrt(static_cast<double>(N))).eval());
    kappas.push_back(1.0);
  }
  const ClosedRing ring(c, h, dropped, kappas, true);
  inst.zeroEigenvalue = dot(dropped.front(), ringApply(h, c, dropped.front()));
  inst.negativeModes = ring.negative();
  Ring mode(static_cast<size_t>(N),
            VectorXd::Ones(f) / std::sqrt(static_cast<double>(f)));
  inst.negativeEigenvalue = lowestRingMode(h, c, mode, 60);
  const long kept = N * f - 1 - m;
  const double logProd =
      static_cast<double>(kept) * std::log(bnh) + 0.5 * ring.logAbsDet();
  inst.logRateTimesZr = -std::log(bnh) +
                        0.5 * std::log(inst.bN / (2.0 * std::numbers::pi *
                                                  inst.betaN * kHbar * kHbar)) -
                        logProd - inst.betaN * inst.ringPotential;

  const auto er = symEigen(hessReactant, false);
  const VectorXd &lr = er.eigenvalues();
  // The m eigenvalues nearest zero are the rigid ones.
  std::vector<long> order(static_cast<size_t>(lr.size()));
  std::iota(order.begin(), order.end(), 0L);
  std::sort(order.begin(), order.end(),
            [&](long a, long b) { return std::abs(lr(a)) < std::abs(lr(b)); });
  std::vector<bool> rigidR(static_cast<size_t>(lr.size()), false);
  for (long k = 0; k < std::min<long>(m, lr.size()); ++k) {
    rigidR[static_cast<size_t>(order[static_cast<size_t>(k)])] = true;
  }
  for (long i = 0; i < lr.size(); ++i) {
    if (!rigidR[static_cast<size_t>(i)] && !(lr(i) > 0.0)) {
      throw std::runtime_error(
          "instantonRate: the reactant Hessian is not positive definite");
    }
  }
  // The rigid modes leave the centroid (k = 0) factor, as the ring's own
  // rigid modes leave its product; both carry them as free particles for
  // k > 0.
  double logZr = -inst.beta * vReactant;
  for (long k = 0; k < N; ++k) {
    const double sk = std::sin(std::numbers::pi * static_cast<double>(k) /
                               static_cast<double>(N));
    for (long i = 0; i < lr.size(); ++i) {
      if (k == 0 && rigidR[static_cast<size_t>(i)]) {
        continue;
      }
      const double l = rigidR[static_cast<size_t>(i)] ? 0.0 : lr(i);
      logZr -= std::log(bnh) + 0.5 * std::log(l + 4.0 * c * sk * sk);
    }
  }
  inst.logZr = logZr;
  inst.logRate = inst.logRateTimesZr - logZr;
  inst.rate = std::exp(inst.logRate) / kTimeUnitSeconds;
  inst.effectiveBarrier =
      -std::log(2.0 * std::numbers::pi * kHbar * inst.beta) / inst.beta -
      inst.logRate / inst.beta;

  if (hessSaddle.size() > 0) {
    const auto ets = symEigen(hessSaddle, false);
    const VectorXd &ls = ets.eigenvalues();
    std::vector<long> orderS(static_cast<size_t>(ls.size()));
    std::iota(orderS.begin(), orderS.end(), 0L);
    std::sort(orderS.begin(), orderS.end(), [&](long a, long b) {
      return std::abs(ls(a)) < std::abs(ls(b));
    });
    std::vector<bool> rigidS(static_cast<size_t>(ls.size()), false);
    for (long k = 0; k < std::min<long>(m, ls.size()); ++k) {
      rigidS[static_cast<size_t>(orderS[static_cast<size_t>(k)])] = true;
    }
    double logRatio = 0.0;
    for (long i = 0; i < lr.size(); ++i) {
      if (!rigidR[static_cast<size_t>(i)]) {
        logRatio += 0.5 * std::log(lr(i));
      }
    }
    // Eigenvalue 0 is the unstable mode, the most negative.
    for (long i = 1; i < ls.size(); ++i) {
      if (!rigidS[static_cast<size_t>(i)]) {
        logRatio -= 0.5 * std::log(std::abs(ls(i)));
      }
    }
    inst.classicalLogRate = logRatio - std::log(2.0 * std::numbers::pi) -
                            inst.beta * (vSaddle - vReactant);
    inst.classicalRate = std::exp(inst.classicalLogRate) / kTimeUnitSeconds;
  }
}

void instantonRate(RateInstanton &inst, const RingBeadHessian &hessian,
                   const MatrixXd &hessReactant, double vReactant,
                   const MatrixXd &hessSaddle, double vSaddle, long rigidModes,
                   long /*denseLimit*/) {
  MatrixXd basis;
  if (rigidModes > 0) {
    const auto er = symEigen(hessReactant);
    const VectorXd &lr = er.eigenvalues();
    std::vector<long> order(static_cast<size_t>(lr.size()));
    std::iota(order.begin(), order.end(), 0L);
    std::sort(order.begin(), order.end(), [&](long a, long b) {
      return std::abs(lr(a)) < std::abs(lr(b));
    });
    const long m = std::min<long>(rigidModes, lr.size());
    basis.resize(lr.size(), m);
    for (long k = 0; k < m; ++k) {
      basis.col(k) = er.eigenvectors().col(order[static_cast<size_t>(k)]);
    }
  }
  instantonRate(inst, hessian, hessReactant, vReactant, hessSaddle, vSaddle,
                basis);
}

double cyclicRingLogAbsDet(double c, const std::vector<MatrixXd> &diag) {
  if (diag.empty()) {
    throw std::invalid_argument("cyclicRingLogAbsDet: no blocks");
  }
  const long f = diag.front().rows();
  std::vector<MatrixXd> h;
  h.reserve(diag.size());
  for (const auto &d : diag) {
    h.push_back(d - 2.0 * c * MatrixXd::Identity(f, f));
  }
  try {
    return ClosedRing(c, h, {}, {}, true).logAbsDet();
  } catch (const std::runtime_error &ex) {
    if (std::string(ex.what()).find("closure") != std::string::npos) {
      return -std::numeric_limits<double>::infinity();
    }
    throw;
  }
}

std::vector<VectorXd> cyclicRingSolve(double c,
                                      const std::vector<MatrixXd> &diag,
                                      const std::vector<VectorXd> &rhs) {
  if (diag.empty() || diag.size() != rhs.size()) {
    throw std::invalid_argument(
        "cyclicRingSolve: one right-hand side per block");
  }
  const long f = diag.front().rows();
  std::vector<MatrixXd> h;
  h.reserve(diag.size());
  for (size_t j = 0; j < diag.size(); ++j) {
    if (rhs[j].size() != f) {
      throw std::invalid_argument("cyclicRingSolve: right-hand side size");
    }
    h.push_back(diag[j] - 2.0 * c * MatrixXd::Identity(f, f));
  }
  const ClosedRing ring(c, h, {}, {}, true);
  if (!std::isfinite(ring.logAbsDet())) {
    throw std::runtime_error("cyclicRingSolve: singular ring");
  }
  return ring.solve(rhs);
}

} // namespace eonc::tunneling
