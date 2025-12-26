# Patient Infographic Review & Update Recommendations

**File:** patient_infographic.png
**Priority:** 8 (Low)
**Status:** Review needed - likely requires minor updates
**Date:** December 26, 2025

---

## Review Checklist

### Elements to Verify

**1. Treatment Recommendations Section**
- [ ] Check if current infographic shows specific drug names
- [ ] If it mentions PI3K/AKT pathway, verify it's up-to-date
- [ ] Check if it references specific medications

**Recommendation:**
- **IF** current infographic shows old/generic treatment info:
  - Update to mention: "Targeted therapy: Alpelisib + Capivasertib"
  - Add: "Clinical trial opportunity available (NCT03602859)"
- **IF** current infographic is generic (no specific drugs mentioned):
  - May not need update, depending on design intent

---

### Elements to ADD (if not present)

**2. Quality Assurance Badge** ⭐ NEW
**Visual:** Small checkmark badge or seal
**Text:** "Quality Checked ✓" or "Advanced QC Validated"
**Location:** Top right corner or near patient ID
**Purpose:** Reassure patient that data quality was rigorously verified

**Example Badge Design:**
```
┌─────────────────┐
│  ✓ QUALITY      │
│  CHECKED        │
│                 │
│  Advanced Data  │
│  Validation     │
└─────────────────┘
```

---

**3. PI3K/AKT/mTOR Pathway Visualization** (if not present)
**Visual:** Simple pathway diagram (patient-friendly version)
**Elements:**
- 3 connected boxes: PI3K → AKT → mTOR
- Red "X" on PTEN (missing brake)
- Drug target icons on PI3K and AKT
- Simple labels: "Cancer growth signals"

**Example Simplified Diagram:**
```
       PTEN (BRAKE) ❌ Missing

           ↓
    ┌──────────────┐
    │    PI3K ↑    │  ← 💊 Alpelisib blocks here
    └──────┬───────┘
           ↓
    ┌──────────────┐
    │    AKT ↑     │  ← 💊 Capivasertib blocks here
    └──────┬───────┘
           ↓
    ┌──────────────┐
    │    mTOR ↑    │
    └──────┬───────┘
           ↓
   Cancer Growth & Resistance
```

---

**4. Treatment Timeline** (if space allows)
**Visual:** Simple horizontal timeline
**Elements:**
- "Analysis Complete" (current)
- "Start Treatment" (next step)
- "Monitor Progress" (ongoing)
- "Response Assessment" (8 weeks)

---

## Recommended Updates by Section

### Section 1: Patient Information
**Current (assumed):**
- Patient ID: PAT001-OVC-2025
- Diagnosis: Stage IV HGSOC
- Status: Platinum-resistant

**Update (add):**
- ✓ Quality Checked badge
- Date of analysis: December 26, 2025

---

### Section 2: Tumor Analysis Findings
**Current (assumed):**
- Genetic mutations listed
- Pathway activation mentioned

**Update (enhance):**
- Add visual: PI3K/AKT/mTOR pathway diagram
- Highlight: PTEN loss (missing brake)
- Visual indicators: ↑ arrows for activated proteins

---

### Section 3: Treatment Recommendations
**Current (unknown - needs review):**
- Generic "targeted therapy" OR
- Specific drug names (may be outdated)

**Update (specific):**
- Primary: Alpelisib + Capivasertib
- Add icons: 💊 for medications
- Add note: "Clinical trial available"
- Visual: Two medication cards side-by-side

**Medication Cards Example:**
```
┌────────────────┐  ┌────────────────┐
│ Alpelisib      │  │ Capivasertib   │
│ (Piqray®)      │  │ (AZD5363)      │
│                │  │                │
│ Blocks PI3K ↑  │  │ Blocks AKT ↑   │
│ FDA Approved   │  │ Phase III      │
└────────────────┘  └────────────────┘
```

---

### Section 4: Next Steps
**Current (assumed):**
- Consult with oncologist
- Discuss treatment options

**Update (specific):**
- Ask about clinical trial NCT03602859
- Meet with nutrition team (blood sugar management)
- Baseline tests needed before starting

---

## Design Recommendations

### Color Scheme
- **Keep patient-friendly:** Soft colors, not clinical
- **Quality badge:** Green (#4CAF50) for trustworthiness
- **Pathway diagram:** Use traffic light colors
  - Green for healthy/normal
  - Red for overactive/problematic
  - Blue for drug targets
- **Medication cards:** Calming blue background

### Typography
- **Patient ID:** Large, bold, easy to read
- **Section headers:** Clear hierarchy
- **Body text:** Minimum 12pt font (readability)
- **Important callouts:** Highlighted boxes or badges

### Icons/Symbols
- ✓ Checkmark: Quality passed, approved
- ↑ Up arrow: Overactive/increased
- ↓ Down arrow: Decreased/inhibited
- ❌ X mark: Missing/deleted
- 💊 Pill icon: Medications
- 📊 Chart icon: Data/analysis
- 🏥 Hospital icon: Clinical trial

---

## Priority Level Assessment

**Current Priority: LOW (8 out of 8)**

**Rationale:**
- Patient infographic is supplementary (nice-to-have, not essential)
- Patient summary HTML provides same information in more detail
- Medication guide provides all treatment info

**When to Prioritize:**
- If patient requests visual summary
- If being used in patient education materials
- If original infographic has outdated/incorrect information

---

## Decision Matrix

### Scenario A: Current infographic is generic
**Action:** Minor update only
- Add quality checkmark badge
- Update date
- Estimated time: 15 minutes

### Scenario B: Current infographic shows specific (old) treatments
**Action:** Moderate update needed
- Replace old drug names with Alpelisib + Capivasertib
- Add clinical trial mention
- Update pathway diagram
- Estimated time: 30-45 minutes

### Scenario C: Current infographic is highly detailed
**Action:** Comprehensive redesign recommended
- Complete update with all new elements
- Add pathway visualization
- Add medication cards
- Estimated time: 1-2 hours

---

## Validation Steps (After Update)

1. **Readability Test:**
   - [ ] All text readable at arm's length (print at full size)
   - [ ] Color contrast sufficient (accessibility)
   - [ ] No medical jargon without explanation

2. **Accuracy Test:**
   - [ ] Drug names spelled correctly
   - [ ] Clinical trial number correct (NCT03602859)
   - [ ] Pathway diagram matches analysis

3. **Consistency Test:**
   - [ ] Matches information in patient_summary.html
   - [ ] Matches medication_guide.html recommendations
   - [ ] Patient ID consistent

4. **Patient-Friendly Test:**
   - [ ] Could a non-medical person understand it?
   - [ ] Are complex concepts simplified appropriately?
   - [ ] Is it reassuring (not alarming) in tone?

---

## Next Steps

**Immediate:**
1. Locate and review current patient_infographic.png file
2. Assess which scenario (A, B, or C) applies
3. Determine if update is necessary based on review

**If Update Needed:**
1. Use specifications from multiomics_resistance_analysis_SPEC.md for technical accuracy
2. Simplify content for patient audience (avoid clinical terminology)
3. Create draft and review with care team before finalizing

**If No Update Needed:**
1. Document review in this file
2. Note: "Current infographic reviewed - no updates required"
3. Move to final regeneration summary

---

## File Format Specifications (If Regenerating)

**Format:** PNG (portable, universally viewable)
**Resolution:** 300 DPI (print quality)
**Dimensions:** 11" × 8.5" (letter size, landscape) or 8.5" × 11" (portrait)
**File Size:** < 5 MB (for easy emailing)
**Color Space:** sRGB (standard for screens)
**Transparency:** No (solid white background)

**Tools for Creation:**
- Adobe Illustrator (professional)
- Canva (user-friendly templates)
- PowerPoint (simple, widely available)
- Inkscape (free, open-source)

---

**Review Status:** ✅ Specification complete
**Update Status:** ⏸️ Pending - awaiting review of current file
**Estimated Time (if update needed):** 15 minutes to 1 hour (depending on scenario)
**Priority:** Low (can be deferred if time-constrained)
