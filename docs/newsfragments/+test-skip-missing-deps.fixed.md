Unit tests SKIP when packaging fixtures or optional engines are missing
(``EON_TEST_SYSTEMS_DIR`` / ``EON_POTENTIALS_PATH`` / nwchemc·cpmdc), so
``meson test`` succeeds for installable client builds without local test data.
