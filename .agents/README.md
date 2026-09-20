# Agent project ledger

This directory is the persistent project memory. Keep it concise enough to inspect
at the beginning of every work session.

- `prompt.md`: scientific specification supplied by the user.
- `state.yaml`: canonical current state and immediate next action.
- `PLAN.md`: staged implementation and production plan.
- `DECISIONS.md`: durable technical and physics decisions.
- `VALIDATION.md`: completed checks, measured resources, and open gates.
- `SECTORS.md`: momentum stars, little groups, and measured block dimensions.

Update the ledger whenever a decision, validation result, scan window, data schema,
or production constraint changes. Store numerical output outside `.agents/` and link
to it from the ledger.
