# Archived Documentation

⚠️ **WARNING:** This directory contains outdated documentation that is no longer maintained.

**Last Updated:** 2026-01-14

---

## Purpose

This archive preserves historical documentation for reference purposes only. Archived documents should **NOT** be used for:
- Current implementation
- Production deployment
- Research analysis
- Clinical decision-making

**For current documentation**, see [docs/INDEX.md](../INDEX.md).

---

## Current Archive Status

**Files Currently Archived:** 0

**Reason:** All documentation was recently reviewed and updated during the 2026-01 documentation reorganization (Phases 1-6). All files have recent "Last Updated" dates (December 2025 or later) and remain actively maintained.

---

## Archiving Criteria

Documents are moved to this archive when they meet **ANY** of these criteria:

### 1. **Outdated (Time-based)**
- Last updated **before October 2025** (6+ months old)
- **AND** not recently validated as still current
- **AND** superseded by newer documentation

### 2. **Marked Obsolete**
- Explicitly marked as "POC only - NOT for production"
- Contains warnings "DEPRECATED" or "DO NOT USE"
- References removed features or discontinued approaches

### 3. **Superseded**
- Content fully replaced by comprehensive newer docs
- Old persona guides replaced by `/docs/for-*/` directories
- Consolidated content (multiple old files → one new file)

### 4. **No Longer Applicable**
- References technologies or approaches no longer in use
- Describes workflows that have been completely redesigned
- Contains instructions for deprecated deployment methods

### Criteria for **NOT** Archiving:
- ❌ Documents with warnings about specific methods (e.g., "not recommended for production" about one deployment option)
- ❌ Work-in-progress notes within active documents
- ❌ Temporary instructions or troubleshooting steps
- ❌ Phase references for active roadmaps ("Phase 1 current, Phase 2 planned")

---

## Archive Structure

```
/docs/archive/
├── README.md (this file)
├── 2025-q3-q4/
│   ├── ARCHIVE_INDEX.md (manifest of archived files)
│   └── [archived files will go here]
└── [future quarter directories]
```

---

## How to Use This Archive

### Before Referencing Archived Content:

1. ✅ **Check current docs first** - See if the topic is covered in active documentation
2. ✅ **Read the disclaimer** - Understand that archived content may be incorrect
3. ✅ **Check archive date** - Understand how old the information is
4. ✅ **Consult ARCHIVE_INDEX.md** - Understand why it was archived

### What You Can Use Archived Content For:

- 📖 **Historical reference** - Understanding how the system evolved
- 📖 **Learning context** - Why certain design decisions were made
- 📖 **Comparison** - How approaches changed over time

### What You Should **NOT** Use It For:

- ❌ **Current implementation** - Use active docs instead
- ❌ **Production deployment** - Archived methods may be insecure
- ❌ **Research citations** - Cite current documentation
- ❌ **Training materials** - Teach current best practices

---

## Archive Contents

### 2025 Q3/Q4 Archive

**Directory:** `2025-q3-q4/`

**Files:** None (all documentation current as of 2026-01-14)

**See:** [2025-q3-q4/ARCHIVE_INDEX.md](2025-q3-q4/ARCHIVE_INDEX.md) for complete manifest

---

## Questions?

**Q: Why was my document archived?**
A: Check [ARCHIVE_INDEX.md](2025-q3-q4/ARCHIVE_INDEX.md) for the specific reason

**Q: I need information from an archived doc**
A: Check if the topic is covered in current docs first. If not, read with caution and verify against current system behavior.

**Q: Can archived docs be restored?**
A: Yes, if they become relevant again. File an issue with justification.

**Q: How often is this archive reviewed?**
A: Quarterly. Extremely old docs (>2 years) may be permanently deleted after notification.

---

## Contact

**Documentation Issues:** https://github.com/lynnlangit/precision-medicine-mcp/issues
**Current Documentation:** [docs/INDEX.md](../INDEX.md)

---

**Related:**
- [Documentation Index](../INDEX.md) - Current active documentation
- [For Developers](../for-developers/README.md) - How to maintain documentation
- [CONTRIBUTING](../../CONTRIBUTING.md) - Documentation contribution guidelines
