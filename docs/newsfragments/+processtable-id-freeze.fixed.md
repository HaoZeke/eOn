Process registration no longer freezes the process-id counter when the
``processtable`` contains duplicate ids. New process ids are
content-addressed with xxHash (``allocate_process_id`` / xxh64 over the
saddle payload and barrier), not ``len(procs)`` or ``max(id)+1``, so
procdata files are not overwritten. Duplicate appends raise.
