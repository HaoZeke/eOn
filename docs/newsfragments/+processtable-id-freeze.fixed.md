Process registration no longer freezes the process-id counter when the
``processtable`` contains duplicate ids. New processes use
``max(id)+1`` (``get_next_process_id``) instead of ``len(procs)``, so
procdata files are not overwritten and silent kinetic collapse is
avoided. Duplicate appends raise; load reloads from disk for allocation.
