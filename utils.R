# --- GSP BENEFICIARY COUNTRIES LIST ---
# Source: Based on U.S. Trade Representative (USTR) designated beneficiary lists.
# This list is comprehensive but static. For a formal publication, one might use
# year-specific lists if available.
cli::cli_h1("Defining GSP Beneficiary Countries")

gsp_country_names <- c(
  "Albania", "Algeria", "Angola", "Argentina", "Armenia", "Azerbaijan", "Bangladesh", 
  "Belize", "Benin", "Bhutan", "Bolivia", "Bosnia and Herzegovina", "Botswana", 
  "Brazil", "Burkina Faso", "Burma", "Burundi", "Cabo Verde", "Cambodia", 
  "Cameroon", "Central African Republic", "Chad", "Comoros", "Congo (Brazzaville)", 
  "Congo (Kinshasa)", "Cote d'Ivoire", "Djibouti", "Ecuador", "Egypt", "Eswatini", 
  "Ethiopia", "Fiji", "Gabon", "Gambia", "Georgia", "Ghana", "Guinea", 
  "Guinea-Bissau", "Guyana", "Haiti", "India", "Indonesia", "Iraq", "Jordan", 
  "Kazakhstan", "Kenya", "Kiribati", "Kosovo", "Kyrgyzstan", "Laos", "Lebanon", 
  "Lesotho", "Liberia", "Madagascar", "Malawi", "Maldives", "Mali", "Mauritania", 
  "Mauritius", "Moldova", "Mongolia", "Montenegro", "Mozambique", "Namibia", "Nepal", 
  "Niger", "Nigeria", "North Macedonia", "Pakistan", "Papua New Guinea", "Paraguay", 
  "Philippines", "Rwanda", "Samoa", "Sao Tome and Principe", "Senegal", "Serbia", 
  "Sierra Leone", "Solomon Islands", "Somalia", "South Africa", "South Sudan", 
  "Sri Lanka", "Sudan", "Suriname", "Tanzania", "Thailand", "Timor-Leste", 
  "Togo", "Tonga", "Tunisia", "Tuvalu", "Uganda", "Ukraine", "Uzbekistan", 
  "Vanuatu", "Yemen", "Zambia", "Zimbabwe"
)

# Convert names to ISO3C codes, which are used in the main dataset
gsp_iso3c_list <- countrycode::countrycode(gsp_country_names, "country.name", "iso3c")
gsp_iso3c_list <- gsp_iso3c_list[!is.na(gsp_iso3c_list)] # Remove any NAs from failed matches

cat("Created list of", length(gsp_iso3c_list), "GSP beneficiary countries (ISO3C).\n")