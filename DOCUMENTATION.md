# Documentation

This document provides an overview of the codebase for generating synthetic MWOs by LLMs and humanising them using rule-based approaches.

## MWO Sentence Generation via LLM

The following functionalities are implemented:

- [`llm_generate.py`](https://github.com/nlp-tlp/Hons24_AllisonLau/blob/main/Generate/llm_generate.py): code to prepare few-shot examples, generate synthetic MWO sentences using GP-4o mini, processing of LM outputs
    - `get_all_paths()`: get all stored paths from json files in `path_patterns` directory
    - `get_generate_prompt()`: prepare prompt for LLM to generate synthetic MWO sentences
    - `get_generate_fewshot()`: prepare few-shot examples for LLM
    - `generate_mwo()`: generate synthetic MWO sentences using LLM (simple)
    - `generate_diverse_mwo()`: generate diverse synthetic MWO sentences using LLM
    - `process_mwo_response()`: process LLM outputs of synthetic MWO sentences
    - `get_samples()`: samples paths from each path type
- [`llm_prompt.py`](https://github.com/nlp-tlp/Hons24_AllisonLau/blob/main/Generate/llm_prompt.py): code to get list of prompt variations, processing of LLM outputs, and paraphrasing the prompts
    - `initialise_prompts()`: get list of prompt variations for LLM
    - `check_similarity()`: check similarity between prompt variations
    - `process_prompt_response()`: process LLM outputs of prompt variations
    - `paraphrase_prompt()`: paraphrase the prompts for LLM
- [`diversity_experiment.ipynb`](https://github.com/nlp-tlp/Hons24_AllisonLau/blob/main/Generate/diversity_experiment.ipynb): experiments for increasing the diversity of the LLM-generated MWO sentences per path
    - Same prompt VS Variations of prompt
    - Single generation VS Batch generation
    - Generation VS Paraphrasing
- You can find the few-shot examples used in the `fewshot_messages` directory.
- Some logs of the LLM-generated MWO sentences can be found in the `mwo_sentences` directory.

## MWO Sentence Humanisation via Rule-based Approach

The following functionalities are implemented:

- [`humanise_experiment.ipynb`](https://github.com/nlp-tlp/Hons24_AllisonLau/blob/main/Humanise/humanise_experiment.ipynb): experiments for humanising synthetic MWO sentences
- [`humanise.py`](https://github.com/nlp-tlp/Hons24_AllisonLau/blob/main/Humanise/humanise.py): rule-based approach for humanising synthetic MWO sentences
    - `initialise_globals()`: initialise global dictionaries for humanisation
    - `load_dictionary()`: load the corrections dictionary for humanisation
    - `shuffle_dictionary()`: shuffle the corrections dictionary
    - `introduce_contractions()`: introduce English contractions to synthetic MWO sentences (50% probability)
    - `introduce_abbreviations()`: introduce abbreviations/jargon to synthetic MWO sentences (40% probability)
    - `rule_introduce_typos()`: introduce up to 3 typos in the synthetic MWO sentences 
    - `humanise_sentence()`: apply the above rules to humanise synthetic MWO sentences


