# Archive Index: 2025 Q3/Q4

**Archive Period:** July 2025 - December 2025
**Archive Date:** January 14, 2026
**Archived By:** Documentation Reorganization (Phases 1-6)

---

## Summary

**Total Files Archived:** 0

**Reason for Empty Archive:**

During the January 2026 documentation reorganization (Phases 1-6), all documentation was comprehensively reviewed:

1. **Phase 1:** Separated test code from test documentation
2. **Phase 2:** Created audience-specific `/docs/for-*/` directories
3. **Phase 3:** Created comprehensive prompt library
4. **Phase 4:** Consolidated overlapping content (server status, costs, HIPAA, risk assessment)
5. **Phase 5:** Reorganized files to logical locations, removed redundant guides
6. **Phase 6:** Assessed all files for archiving

**Findings:**
- ✅ All documentation has "Last Updated" dates from December 2025 or January 2026
- ✅ No files explicitly marked as "deprecated" or "obsolete"
- ✅ All Phase references are for active roadmaps or completed milestones
- ✅ Security warnings refer to specific deployment methods, not entire documents
- ✅ All content actively maintained and relevant

**Conclusion:** No files met the archiving criteria. All documentation remains current and valid.

---

## Files Reviewed for Archiving

The following files were explicitly checked against archiving criteria:

| File | Last Updated | Status | Decision |
|------|--------------|--------|----------|
| `/docs/deployment/security.md` | Dec 29, 2025 | Active | ✅ Keep - Recent, comprehensive security guide |
| `/docs/for-hospitals/compliance/hipaa.md` | Dec 2025 | Active | ✅ Keep - Current HIPAA compliance documentation |
| `/docs/architecture/servers.md` | Jan 10, 2026 | Active | ✅ Keep - Consolidated server status (Phase 4) |
| `/docs/for-*/README.md` (6 files) | Jan 14, 2026 | Active | ✅ Keep - New comprehensive audience guides (Phase 2) |
| `/docs/for-developers/automation-guides/for-*.md` (6 files) | Various | **Deleted** | 🗑️ Removed in Phase 5 (superseded by for-*/ directories) |

### Files Deleted (Not Archived)

**Reason:** These files were removed during Phase 5 because they were **superseded by better content**, not because they were outdated. The information was migrated to comprehensive directories.

| Deleted File | Superseded By | Migration Date |
|--------------|---------------|----------------|
| `docs/for-developers/automation-guides/for-developers.md` | `docs/for-developers/README.md` | Jan 14, 2026 |
| `docs/for-developers/automation-guides/for-researchers.md` | `docs/for-researchers/README.md` | Jan 14, 2026 |
| `docs/for-developers/automation-guides/for-clinicians.md` | `docs/for-hospitals/citl-workflows/CITL_WORKFLOW_GUIDE.md` | Jan 14, 2026 |
| `docs/for-developers/automation-guides/for-bioinformaticians.md` | `docs/for-researchers/README.md` | Jan 14, 2026 |
| `docs/for-developers/automation-guides/for-patients.md` | `docs/for-patients/README.md` | Jan 14, 2026 |
| `do../../for-developers/ADD_NEW_MODALITY_SERVER.md` | `docs/for-developers/ADD_NEW_MODALITY_SERVER.md` | Jan 14, 2026 |

**Note:** These files were **deleted, not archived** because:
- Content was fully migrated to new locations
- New versions are more comprehensive
- Keeping old versions would create confusion
- No historical value (recent content, not old approaches)

---

## Archiving Criteria Applied

### Criteria 1: Time-based Outdated
**Assessment:** ❌ No files met this criteria
- All docs last updated December 2025 or later
- All content recently validated during reorganization

### Criteria 2: Marked Obsolete
**Assessment:** ❌ No files met this criteria
- No files explicitly marked "DEPRECATED" or "DO NOT USE"
- Warnings found were about specific methods within documents, not entire docs

### Criteria 3: Superseded
**Assessment:** ✅ 6 files met this criteria, but were **deleted** not archived
- Reason: Recent content with no historical value
- Content fully migrated to new comprehensive directories
- Archiving would imply they had historical value worth preserving

### Criteria 4: No Longer Applicable
**Assessment:** ❌ No files met this criteria
- All deployment methods still supported
- All technologies still in use
- All workflows still applicable

---

## Future Archiving Guidelines

### When to Archive vs. Delete

**Archive if:**
- Document represents significant historical approach
- Contains valuable context about design decisions
- References abandoned but interesting experiments
- Last updated >6 months ago AND superseded

**Delete if:**
- Recent content fully migrated elsewhere
- Duplicate/redundant with no unique information
- Pure administrative docs (old meeting notes, etc.)
- No historical or reference value

### Recommended Archive Schedule

**Quarterly Review** (every 3 months):
1. Check all "Last Updated" dates
2. Identify docs not updated in 6+ months
3. Verify if content is still accurate
4. Archive or update as needed

**Next Review:** April 2026

---

## Restoration Process

If an archived document becomes relevant again:

1. Open GitHub issue explaining why restoration is needed
2. Review document for accuracy
3. Update content to current standards
4. Move back to appropriate `/docs/` location
5. Update "Last Updated" date
6. Remove from archive index

---

## Contact

**Questions about this archive:** File an issue at https://github.com/lynnlangit/precision-medicine-mcp/issues

**Current documentation:** See [/docs/INDEX.md](../../INDEX.md)

---

**Last Updated:** 2026-01-14
**Next Review:** 2026-04 (Q2)
