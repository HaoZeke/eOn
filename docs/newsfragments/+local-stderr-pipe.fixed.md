Local communicator redirects client stderr to `stderr.dat` beside `stdout.dat`
instead of an undrained `PIPE`. Clients that wrote more than one OS pipe buffer
(~64 KiB) previously blocked forever in `anon_pipe_write` while the parent only
polled `p.poll()`, so large aKMC saddle searches never returned results.
