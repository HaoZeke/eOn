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
#pragma once

#include "../../Potential.h"
#include <memory>
#include <string>
#include <vector>

// --- Constants for i-PI protocol and unit conversions ---
// (Must match the values used in SocketIPIPot.cpp)

// Message length for i-PI header
constexpr int IPI_MSG_LEN = 12;
// Bohr in Angstrom
constexpr double BOHR_IN_ANGSTROM = 0.52917721092;
// Hartree in eV
constexpr double HARTREE_IN_EV = 27.21138602;

namespace {
// Helper to get element symbols from atomic numbers
// (Must be in the header if it's not in the .cpp)
const char *elementArray[] = {
    "Unknown", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na",      "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti",      "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As",      "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru",      "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe",
    "Cs",      "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy",      "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os",
    "Ir",      "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra",      "Ac", "Th", "Pa", "U",  nullptr};
char const *atomicNumber2symbol(int n) { return elementArray[n]; }
} // namespace

/**
 * @brief A potential that acts as a CLIENT to an i-PI server.
 *
 * Implements the i-PI socket protocol to communicate with an i-PI executable
 * that has been launched in server mode (e.g., using <ffsocket mode='driver'>).
 * This class manages a persistent connection TO the i-PI server.
 * Supports both TCP/IP and UNIX domain sockets.
 */
class IPIPot : public Potential {
public:
  /**
   * @brief Constructor. Reads config, creates socket, and connects to i-PI.
   */
  explicit IPIPot(std::shared_ptr<Parameters> p);

  /**
   * @brief Destructor. Sends "EXIT" to i-PI and closes the socket.
   */
  ~IPIPot() override;

  /**
   * @brief The method called by Eon to compute forces and energy.
   *
   * This method implements the full i-PI client-side handshake:
   * 1. Checks server status, handling NEEDINIT if required.
   * 2. Sends atomic positions and cell (POSDATA).
   * 3. Polls until the calculation is complete (HAVEDATA).
   * 4. Retrieves the mean force, mean energy, and variance (GETFORCE).
   *
   * @param N The number of atoms.
   * @param R Pointer to the position array [x1, y1, z1, ...].
   * @param atomicNrs Pointer to the array of atomic numbers.
   * @param F Pointer to the force array to be populated (mean force).
   * @param U Pointer to the potential energy value to be set (mean energy).
   * @param variance Pointer to the energy variance (from i-PI's sampling).
   * @param box Pointer to the 3x3 simulation box matrix.
   */
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  // --- Private Methods ---

  /**
   * @brief Creates the socket and connects to the i-PI server.
   */
  void connect_to_server();

  /**
   * @brief Sends a fixed-length (IPI_MSG_LEN) header to the server.
   */
  void send_header(const char *msg);

  /**
   * @brief Receives a fixed-length (IPI_MSG_LEN) header from the server.
   */
  void recv_header(char *buffer);

  /**
   * @brief Sends an exact number of bytes (n_bytes) to the server.
   */
  void send_exact(const void *buffer, size_t n_bytes);

  /**
   * @brief Receives an exact number of bytes (n_bytes) from the server.
   */
  void recv_exact(void *buffer, size_t n_bytes);

  // --- Member Variables ---
  int sock_fd; // Socket file descriptor (for the client)
  bool is_connected;
  bool unix_socket_mode;
  std::string server_address; // Host for TCP, path for UNIX
  int port;

  // Cache for atom symbols, needed for the one-time INIT message
  std::vector<std::string> atom_symbols_cache;
};
