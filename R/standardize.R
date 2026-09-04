#' Standardize FAERS Quarterly Data for Preferred Term and drug names
#' @param object A [FAERSascii] object.
#' @param ... Other arguments passed to specific methods.
#' @return A [FAERSascii] object.
#' @seealso [MedDRA]
#' @examples
#' #' # you must change `dir`, as the files included in the package are sampled
#' data <- faers(c(2004, 2017), c("q1", "q2"),
#'     dir = system.file("extdata", package = "faers"),
#'     compress_dir = tempdir()
#' )
#' \dontrun{
#' # you should replace `meddra_path` with yours
#' data <- faers_standardize(data, meddra_path)
#' }
#' @export
#' @name faers_standardize
methods::setGeneric("faers_standardize", function(object, ...) {
    standardGeneric("faers_standardize")
})

#' @param meddra_path A string, define the path of MedDRA directory.
#' @inheritParams meddra
#' @export
#' @method faers_standardize FAERSascii
#' @rdname faers_standardize
methods::setMethod("faers_standardize", "FAERSascii", function(object, meddra_path, add_smq = FALSE) {
    if (!is.null(object@db)) {
        return(standardize_db(object, meddra_path, add_smq = add_smq))
    }
    # standardize PT terms
    # for indi
    assert_string(meddra_path)
    meddra_data <- meddra(meddra_path, add_smq = add_smq, primary_soc = TRUE)
    # https://stackoverflow.com/questions/70181149/is-a-saved-and-loaded-data-table-with-qs-a-correct-data-table
    # fix error: when load a saved FAERS object, don't change by reference
    cli::cli_alert("standardize {.field Preferred Term} in indi")
    object@data$indi <- data.table::as.data.table(object@data$indi)

    meddra_cols <- c("meddra_hierarchy_idx", "meddra_hierarchy_from", "meddra_code", "meddra_pt")
    object@data$indi[, (meddra_cols) := meddra_standardize_pt(
        clean_indi_pt(indi_pt, meddra_data@hierarchy), # nolint
        meddra_data@hierarchy
    )]

    cli::cli_alert("standardize {.field Preferred Term} in reac")
    object@data$reac <- data.table::as.data.table(object@data$reac)
    object@data$reac[, (meddra_cols) := meddra_standardize_pt(
        clean_reac_pt(pt, meddra_data@hierarchy), # nolint
        meddra_data@hierarchy
    )]
    object@meddra <- meddra_data
    object@standardization <- TRUE
    object
})

clean_indi_pt <- function(x, hierarchy) {
    x <- str_replace_all(str_trim(toupper(x)), "\\s+", " ")
    code <- data.table::fcase(
        x == "ACID REFLUX", "10017885",
        x == "ACUTE ISCHEMIC STROKE", "10061256",
        x == "AEROMONA INFECTION", "10054205",
        x == "ARTRIAL FIBRILLATION", "10003658",
        x == "ATRIAL FILBRILLATION", "10003658",
        x == "AUTIOIMMUNE INDUCED RASH", "10075689", # Autoimmune dermatitis
        x == "B LYMPHOBLASTIC LEUKEMIA", "10054448",
        x == "BACTERIAL PNEUMONIA", "10060946",
        x == "BIPOLAR DISORDER II", "10004940",
        x == "BLODD PRESSURE", "10005727",
        x == "CARDIAC CATH", "10007527",
        x == "CARDIC DISORDER", "10061024",
        x == "CORONARY ARTERY DISEASE/HYPERTENSION", "10020772", #  Hypertension
        x == "CHEMOTHERAPY/RECTOSIGMOID CANCER", "10038093",
        x == "CHRONIC NERVE PAIN", "10029181",
        x == "COMPLEX PARTIAL EPILEPSY", "10010145",
        x == "COROARY ARTERY STENT PLACEMENT", "10052086",
        x == "CORONARY STENT INSERTION", "10052086",
        x == "CROHNS DISEASE", "10011401",
        x == "DEPRESIION", "10012378",
        x == "DEPRESSION NEC", "10012378",
        x == "DEPRESSON", "10012378",
        x == "DISBACTERIOSIS", "10064389",
        x == "DRUG", "10063370", # 10063370 Drug therapy
        x == "DRUG INDUCED LIVER INJURY", "10072268",
        x == "DRUG KNOWN FOR UNKNOWN INDICATION", "10057097",
        x == "DRUG USE", "10063370", # 10063370 Drug therapy
        x == "DRUG USE FO RUNKNOWN INDICATION", "10057097",
        x == "DRUG USE UNKNOWN INDICATION", "10057097",
        x == "DRUG USED FOR UNKNOWN INDICATION", "10057097",
        x == "EVAN'S SYNDROME", "10053873",
        x == "EXUDATIVE SENILE MACULAR DEGENERATIVE OF RETINA", "10015902",
        x == "FOLLICLE-STIMULATING HORMONE DEFICIENCY", "10071084",
        x == "GASTRO-JEJUNOSTOMY", "10017882",
        x == "GBM", "10018336",
        x == "GENERAL ANXIETY DISORDER", "10018075",
        x == "H1NI INFLUENZA", "10069767",
        x == "HEADACHE/PROPHYLAXIS", "10019211",
        x == "HEADACHES", "10019211",
        x == "HEAVY BLEEDING", "10005103",
        x == "HELICOBACTER PYLORI PROPHYLAXIS", "10054263",
        x == "HEPATOBILLIARY DISORDER PROPHYLAXIS", "10081385",
        x == "HER-2 POSITIVE BREAST CANCER", "10065430",
        x == "HER-2 POSITIVE GASTRIC CANCER", "10066896",
        x == "HER-2 PROTEIN OVEREXPRESSION", "10075638",
        x == "HUNTINGTONS DISEASE", "10070668",
        x == "HYPOTHALAMO-PITUITARY DISORDERS", "10021111",
        x == "HYPOTHYROID", "10021114",
        x == "INFECTIVE ENDOCARDITIS", "10014678",
        x == "INFECTED MOLE", "10027806", # Mole of skin
        x == "INFLAMATION", "10061218", # 10061218 Inflammation
        x == "INHALATION", "10052996", # 10052996 Inhalation therapy
        x == "INR INCREAESD", "10022402",
        x == "INTERVERTEBRAL DISKITIS", "10060738",
        x == "INTRACTABLE SPASTICITY", "10041416", # Spasticity
        x == "IRREGULAR HEARTBEAT", "10019323",
        x == "LENNOX--GASTAUT SYNDROME", "10048816",
        x == "LUMBAR SPONDYLOSIS WITH SCOLIOSIS AND ARTHRITIS", "10025007",
        x == "MENINGEOMAS SURGERY", "10053765",
        x == "METASTATIC CHOLANGIOCARCINOMA", "10077846",
        x == "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS INFECTION", "10027508",
        x == "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS SEPSIS", "10058867",
        x == "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS TEST", "10053429",
        x == "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS TEST POSITIVE", "10053427",
        x == "NECROTIZING PNEUMONIA", "10049271",
        x == "OSEOPOROSIS", "10031282",
        x == "OSTEOMYLITIS", "10031252",
        x == "OSTEOPORISIS", "10031282",
        x == "PAH", "10064911",
        x == "PANIC DISORDER/ DEPRESSION", "10033666",
        x == "PRECAUTIONARY MEASURE", "10036898", # 10036898 Prophylaxis
        x == "PRODUCT USED FOR UNKNOWN INDCATION", "10070592",
        x == "PTSD", "10036316",
        x == "PULMONARY ARTERY HYPERTENSION", "10064911",
        x == "REACTIVE AIRWAY DISEASE", "10037993",
        x == "REGULATE HEART RATE", "10019304", # 10019304 Heart rate irregular
        x == "SCHIZOPRENIA", "10039626",
        x == "SEIZURE DISORDER", "10039906",
        x == "SIEZURE", "10039906",
        x == "SKIN AND SOFT TISSUE INFECTION", "10081417",
        x == "SMOKING CESSATION", "10008374",
        x == "SPINAL LUMBAR DISORDER", "10061368",
        x == "STAGE IV NON-SMALL CELL", "10029522",
        x == "STAGE IV RECTAL ADENOCARCINOMA", "10038029",
        x == "STREPTOCOCCAL IDENTIFICATION TEST", "10067006",
        x == "STROKE PREVENTION", "10081388", # Nervous system disorder prophylaxis
        x == "SUNBURN PROPHYLAXIS", "10081391", # Skin disorder prophylaxis
        x == "SUPERIOR VENA CAVAL OCCLUSION", "10042568",
        x == "THYROID CONDITION", "10043710",
        x == "THYROID HORMONE REPLACEMENT", "10068076",
        x == "TOTAL KNEE ARTHROPLASTY", "10003398",
        x == "TRAVELLER'S DIARRHEA", "10044552",
        x == "TYPE II DIABETES", "10045242"
        # following items were not mapped
        # "LLT" : drug  ACTONEL (RISEDRONATE SODIUM) TABLET, 150MG ?
    )
    operated_idx <- !is.na(code)
    x[operated_idx] <- meddra_map_code_into_names(
        codes = code[operated_idx], hierarchy
    )
    x
}

clean_reac_pt <- function(x, hierarchy) {
    x <- str_replace_all(str_trim(toupper(x)), "\\s+", " ")
    code <- data.table::fcase(
        x == "ANO-RECTAL STENOSIS", "10002581",
        x == "HER-2 POSITIVE BREAST CANCER", "10065430",
        x == "STAPHYLOCOCCAL IDENTIFICATION TEST POSITIVE", "10067140",
        x == "STREPTOCOCCAL IDENTIFICATION TEST", "10067006",
        x == "STREPTOCOCCAL SEROLOGY", "10059987",
        x == "ANO-RECTAL ULCER", "10002582",
        x == "BLASTIC PLASMACYTOID DENDRITRIC CELL NEOPLASIA", "10075460",
        x == "CORNELIA DE-LANGE SYNDROME", "10077707",
        x == "FRONTAL SINUS OPERATIONS", "10017379",
        x == "HER-2 POSITIVE GASTRIC CANCER", "10066896",
        x == "MAXILLARY ANTRUM OPERATIONS", "10026950",
        x == "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS TEST NEGATIVE", "10053428",
        x == "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS TEST POSITIVE", "10053427",
        x == "PAROVARIAN CYST", "10052456",
        x == "STAPHYLOCOCCAL IDENTIFICATION TEST NEGATIVE", "10067005",
        x == "STREPTOCOCCAL IDENTIFICATION TEST POSITIVE", "10067004",
        x == "AEROMONA INFECTION", "10054205",
        x == "DISBACTERIOSIS", "10064389",
        x == "GASTRO-ENTEROSTOMY", "10017873",
        x == "HER-2 PROTEIN OVEREXPRESSION", "10075638",
        x == "HYPOTHALAMO-PITUITARY DISORDERS", "10021111",
        x == "SUPERIOR VENA CAVAL OCCLUSION", "10058988",
        x == "AEROMONA INFECTION", "10054205",
        x == "CAPNOCYTOPHAGIA INFECTION", "10061738",
        x == "DISBACTERIOSIS", "10064389",
        x == "EAGLES SYNDROME", "10066835",
        x == "EVAN'S SYNDROME", "10053873",
        x == "GASTRO-INTESTINAL FISTULA", "10071258",
        x == "GLYCOPEPTIDE ANTIBIOTIC RESISTANT STAPHYLOCOCCAL AUREUS INFECTION", "10052101",
        x == "HEPATOBILLIARY DISORDER PROPHYLAXIS", "10081385",
        x == "MENINGEOMAS SURGERY", "10053765",
        x == "METHICILLIN-RESISTANT STAPHYLOCOCCAL AUREUS TEST", "10053429",
        x == "SPHENOID SINUS OPERATIONS", "10041508",
        x == "STREPTOCOCCAL SEROLOGY POSITIVE", "10059988",
        x == "SUPERIOR VENA CAVAL STENOSIS", "10064771",
        x == "GASTRO-JEJUNOSTOMY", "10017882",

        # self-added
        x == "IMMUNE-MEDIATED ADRENAL INSUFICIENCY", "10085547",
        x == "ORAL APPLIANCE", "10085270",
        x == "VAGINAL RING", "10082353",
        x == "VAGINAL CUFF", "10088846",
        x == "MOREL-LAVELLEE SEROMA", "10088873",
        x == "RABSON MENDENHALL SYNDROME", "10088742"
    )
    operated_idx <- !is.na(code)
    x[operated_idx] <- meddra_map_code_into_names(
        codes = code[operated_idx], hierarchy
    )
    x
}

########################################################
faers_standardize_drug <- function(terms, force = FALSE, url = NULL) {
    athena_standardize_drug(terms = terms, force = force, url = url)
}

utils::globalVariables(c("cleaned_pt", "indi_pt", "pt"))
