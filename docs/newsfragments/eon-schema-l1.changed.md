Move full job-config pydantic models into ``eon-schema`` (``eon_schema.config``).
``eon.schema`` and ``pyeonclient.models`` re-export the shared package so both
eon-akmc and pyeonclient share one schema surface.
