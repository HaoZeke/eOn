#include "IPIPot.h"

#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>


using std::this_thread::sleep_for;

/**
 * @brief Constructor for IPIPot.
 *
 * Reads i-PI socket parameters, creates the socket, and connects to the
 * i-PI server. Performs the initial handshake.
 */
IPIPot::IPIPot(std::shared_ptr<Parameters> p)
    : Potential(PotType::IPI, p),
      sock_fd(-1),
      is_connected(false) {

  // --- IMPORTANT ---
  // You must add a 'socket_ipi_options' struct to your Parameters.h
  // similar to 'socket_nwchem_options'
  //
  // struct SocketIPiOptions {
  //     bool unix_socket_mode = true;
  //     std::string unix_socket_path = "/tmp/ipi_eon_socket";
  //     std::string host = "127.0.0.1";
  //     int port = 31415;
  // };
  //
  // --- END IMPORTANT ---

  unix_socket_mode = p->socket_ipi_options.unix_socket_mode;

  if (unix_socket_mode) {
    server_address = p->socket_ipi_options.unix_socket_path;
    port = -1;
    std::cout << "IPIPot: Initializing in UNIX mode." << std::endl;
    std::cout << "Connecting to i-PI server at: " << server_address
              << std::endl;
  } else {
    server_address = p->socket_ipi_options.host;
    port = p->socket_ipi_options.port;
    std::cout << "IPIPot: Initializing in TCP mode." << std::endl;
    std::cout << "Connecting to i-PI server at: " << server_address << ":"
              << port << std::endl;
  }

  // Connect to the i-PI server
  connect_to_server();
  is_connected = true;
  std::cout << "Successfully connected to i-PI server." << std::endl;

  // Perform initial handshake
  char status_buffer[IPI_MSG_LEN + 1] = {0};
  send_header("STATUS");
  recv_header(status_buffer);
  if (std::string(status_buffer) != "READY" &&
      std::string(status_buffer) != "NEEDINIT") {
    throw std::runtime_error("i-PI server not ready after connection. Sent: " +
                             std::string(status_buffer));
  }
  std::cout << "i-PI server is " << status_buffer << std::endl;
}

/**
 * @brief Destructor.
 *
 * Sends an "EXIT" message to the i-PI server and closes the socket.
 */
IPIPot::~IPIPot() {
  if (is_connected) {
    std::cout << "Closing connection to i-PI server..." << std::endl;
    try {
      send_header("EXIT");
    } catch (...) {
      // Ignore errors during shutdown
    }
  }
  if (sock_fd >= 0)
    ::close(sock_fd);
}

/**
 * @brief Main force call implementation (see header for details).
 */
void IPIPot::force(long N, const double *R, const int *atomicNrs,
                         double *F, double *U, double *variance,
                         const double *box) {
  if (!is_connected) {
    throw std::runtime_error("i-PI server is not connected.");
  }

  // Check server status
  char status_buffer[IPI_MSG_LEN + 1] = {0};
  send_header("STATUS");
  recv_header(status_buffer);

  // Handle INIT if i-PI requests it
  if (std::string(status_buffer) == "NEEDINIT") {
    send_header("INIT");

    // Cache atom symbols if this is the first call
    if (atom_symbols_cache.empty()) {
      atom_symbols_cache.reserve(N);
      for (long i = 0; i < N; ++i) {
        atom_symbols_cache.emplace_back(atomicNumber2symbol(atomicNrs[i]));
      }
    }

    // Build the extra string: "ATOMS: H O H ..."
    std::string extra_string = "ATOMS: ";
    for (long i = 0; i < N; ++i) {
      extra_string += " " + atom_symbols_cache[i];
    }
    // Add comment for clarity
    extra_string += " # EON Client";

    int32_t init_payload[2];
    init_payload[0] = 0; // bead_index (we only send one bead)
    init_payload[1] = extra_string.length() + 1; // +1 for null terminator

    send_exact(&init_payload, sizeof(init_payload));
    send_exact(extra_string.c_str(), init_payload[1]); // Send with null term

    // Re-check status
    send_header("STATUS");
    recv_header(status_buffer);
  }

  if (std::string(status_buffer) != "READY") {
    throw std::runtime_error("i-PI server not ready for new positions! Sent: " +
                             std::string(status_buffer));
  }

  // --- Send Position Data ---
  // i-PI expects Angstrom, so no conversion of R is needed.
  // We send identity matrices for cell/inv_cell for simplicity,
  // assuming non-periodic calculations.
  // If you need PBCs, you must calculate and send the *real*
  // transposed cell and inverse cell here.
  double invcell_T[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double cell_T[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  // Or, if `box` is the 9-element cell matrix (row-major):
  // cell_T[0]=box[0]; cell_T[1]=box[3]; cell_T[2]=box[6]; // Col 1
  // cell_T[3]=box[1]; cell_T[4]=box[4]; cell_T[5]=box[7]; // Col 2
  // cell_T[6]=box[2]; cell_T[7]=box[5]; cell_T[8]=box[8]; // Col 3
  // ... and you would need to compute invcell_T ...

  send_header("POSDATA");
  int32_t nat = N;
  send_exact(cell_T, sizeof(cell_T));
  send_exact(invcell_T, sizeof(invcell_T));
  send_exact(&nat, sizeof(nat));
  send_exact(R, N * 3 * sizeof(double)); // Send positions in Angstrom

  // --- Poll for results ---
  while (true) {
    send_header("STATUS");
    recv_header(status_buffer);
    if (std::string(status_buffer) == "HAVEDATA") {
      break;
    }
    // A small sleep to prevent busy-waiting
    sleep_for(10ms);
  }

  // --- Request and receive results ---
  send_header("GETFORCE");
  recv_header(status_buffer);
  if (std::string(status_buffer) != "FORCEREADY") {
    throw std::runtime_error("Expected FORCEREADY, got " +
                             std::string(status_buffer));
  }

  // Unpack the results payload
  double energy_ha;
  int32_t nat_back;
  std::vector<double> forces_ha_bohr(N * 3);
  double virial_ha[9];
  int32_t extra_len;

  recv_exact(&energy_ha, sizeof(energy_ha));
  recv_exact(&nat_back, sizeof(nat_back));
  if (nat_back != N)
    throw std::runtime_error("Atom count mismatch from i-PI");
  recv_exact(forces_ha_bohr.data(), forces_ha_bohr.size() * sizeof(double));
  recv_exact(&virial_ha, sizeof(virial_ha));
  recv_exact(&extra_len, sizeof(extra_len));

  std::vector<char> extra_buf;
  if (extra_len > 0) {
    extra_buf.resize(extra_len);
    recv_exact(extra_buf.data(), extra_len);
  }

  // --- Convert results to EON units ---
  // i-PI sends mean energy in Hartree
  *U = energy_ha * HARTREE_IN_EV;

  // i-PI sends mean forces in Hartree/Bohr
  for (int i = 0; i < N * 3; ++i) {
    F[i] = forces_ha_bohr[i] * (HARTREE_IN_EV / BOHR_IN_ANGSTROM);
  }

  // Parse variance from the extra string
  // We assume the *first* value in the string is the POTENTIAL ENERGY
  // VARIANCE in (Hartree^2).
  if (extra_len > 0) {
    std::string extra_str(extra_buf.data());
    std::stringstream ss(extra_str);
    double variance_ha2;
    if (ss >> variance_ha2) {
      // Convert from (Hartree^2) to (eV^2)
      *variance = variance_ha2 * (HARTREE_IN_EV * HARTREE_IN_EV);
    } else {
      *variance = 0.0; // Failed to parse
    }
  } else {
    *variance = 0.0;
  }
}

// =================================================
// Private Helper Methods for Socket Communication
// =================================================

/**
 * @brief Creates a socket and connects to the i-PI server.
 */
void IPIPot::connect_to_server() {
  int domain = unix_socket_mode ? AF_UNIX : AF_INET;
  sock_fd = socket(domain, SOCK_STREAM, 0);
  if (sock_fd < 0) {
    throw std::runtime_error("Failed to create socket.");
  }

  if (unix_socket_mode) {
    sockaddr_un sock_addr{};
    sock_addr.sun_family = AF_UNIX;
    strncpy(sock_addr.sun_path, server_address.c_str(),
            sizeof(sock_addr.sun_path) - 1);

    socklen_t addr_len =
        sizeof(sock_addr.sun_family) + strlen(sock_addr.sun_path);
    if (::connect(sock_fd, (struct sockaddr *)&sock_addr, addr_len) < 0) {
      perror("UNIX connect failed");
      throw std::runtime_error("Failed to connect to i-PI UNIX socket.");
    }
  } else { // TCP Mode
    sockaddr_in sock_addr{};
    sock_addr.sin_family = AF_INET;
    sock_addr.sin_addr.s_addr = inet_addr(server_address.c_str());
    sock_addr.sin_port = htons(port);

    if (::connect(sock_fd, (struct sockaddr *)&sock_addr, sizeof(sock_addr)) <
        0) {
      perror("TCP connect failed");
      throw std::runtime_error("Failed to connect to i-PI TCP socket.");
    }
  }
}

/**
 * @brief Sends a fixed-length header message.
 */
void IPIPot::send_header(const char *msg) {
  char buffer[IPI_MSG_LEN] = {0};
  strncpy(buffer, msg, IPI_MSG_LEN);
  send_exact(buffer, IPI_MSG_LEN);
}

/**
 * @brief Receives a fixed-length header message.
 */
void IPIPot::recv_header(char *buffer) {
  recv_exact(buffer, IPI_MSG_LEN);
  buffer[IPI_MSG_LEN] = '\0'; // Null-terminate
  // Trim trailing whitespace
  for (int i = IPI_MSG_LEN - 1; i >= 0 && isspace((unsigned char)buffer[i]);
       --i) {
    buffer[i] = '\0';
  }
}

/**
 * @brief Sends an exact number of bytes.
 */
void IPIPot::send_exact(const void *buffer, size_t n_bytes) {
  size_t sent = 0;
  while (sent < n_bytes) {
    ssize_t n = ::send(sock_fd, (const char *)buffer + sent, n_bytes - sent, 0);
    if (n <= 0) {
      throw std::runtime_error(
          "send_exact failed: connection closed or error.");
    }
    sent += n;
  }
}

/**
 * @brief Receives an exact number of bytes.
 */
void IPIPot::recv_exact(void *buffer, size_t n_bytes) {
  size_t recvd = 0;
  while (recvd < n_bytes) {
    ssize_t n = ::recv(sock_fd, (char *)buffer + recvd, n_bytes - recvd, 0);
    if (n <= 0) {
      throw std::runtime_error(
          "recv_exact failed: connection closed or error.");
    }
    recvd += n;
  }
}
