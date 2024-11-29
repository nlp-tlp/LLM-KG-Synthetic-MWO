import os
import re
import csv
import random
import nltk
from openai import OpenAI
from dotenv import load_dotenv
from nltk.corpus import cmudict
from Levenshtein import distance as levenshtein_distance

# Global variables
CONTRACTIONS_DICT = {}      # {expand: [contractions]}
ABBREVIATIONS_DICT = {}     # {original: [variations]}
KEYBOARD_DICT = {}          # {key: [adjacent]}
CMU_DICT = {}               # CMU pronouncing dictionary

def initialise_globals(dirpath):
    global CONTRACTIONS_DICT, ABBREVIATIONS_DICT, KEYBOARD_DICT, CMU_DICT
    path = os.path.join(dirpath, 'data', 'Corrections')
    
    # Load contractions {expand: [contractions]}
    CONTRACTIONS_DICT = load_dictionary(os.path.join(path, 'contractions.csv'))

    # Load abbreviations {original: [variations]}
    ABBREVIATIONS_DICT = load_dictionary(os.path.join(path, 'abbreviations.csv'))

    # Load keyboard adjacent letters {key: [adjacent]}
    KEYBOARD_DICT = load_dictionary(os.path.join(path, 'keyboard.csv'))
    
    # Load CMU pronouncing dictionary
    nltk.download('cmudict')
    CMU_DICT = cmudict.dict()

# Load abbreviations dictionary
def load_dictionary(file):
    """" Load a dictionary from a CSV file where the 
        first column is the original word and the 
        second column is the variation. """
    dictionary = {}
    with open(file, 'r') as f:
        reader = csv.reader(f)
        next(reader) # Ignore header
        for row in reader:
            original = row[0]
            variations = row[1]
            if original in dictionary:
                dictionary[original].append(variations)
            else:
                dictionary[original] = [variations]
    dictionary = dict(sorted(dictionary.items()))
    return dictionary

# Shuffle the dictionary
def shuffle_dictionary(d):
    """ Shuffle the dictionary """
    items = list(d.items())  # Convert dictionary to list of items
    random.shuffle(items)    # Shuffle the list
    return dict(items)       # Convert list back to dictionary

# Introduce contractions in a sentence (default probability=0.5)
def introduce_contractions(sentence, chance=0.5):
    """ Introduce contractions in a sentence """
    contractions = shuffle_dictionary(CONTRACTIONS_DICT) # Shuffle
    for expanded, contracted in contractions.items():
        pattern = r'\b' + expanded + r'\b' # Match whole word
        if re.search(expanded, sentence, re.I) and random.random() < chance: # Case-insensitive
            sentence = re.sub(pattern, random.choice(contracted), sentence, flags=re.I)
    return sentence

# Introduce abbreviations in a sentence (default probability=0.4)
def introduce_abbreviations(sentence, chance=0.4):
    """ Introduce abbreviations in a sentence """
    abbreviations = shuffle_dictionary(ABBREVIATIONS_DICT) # Shuffle
    for original, variations in abbreviations.items():
        pattern = r'\b' + original + r'\b' # Match whole word
        # Check if original word is in sentence
        if re.search(original, sentence, re.I) and random.random() < chance: # Case-insensitive
            variation = random.choice(variations)
            variation = add_periods(original, variation)
            sentence = re.sub(pattern, variation, sentence, flags=re.I)
    return sentence

# Add periods to abbreviations if condition is met (default probability=0.08)
def add_periods(original_word, abbreviated_word, chance=0.08):
    """ Add periods to an abbreviation if the condition is met """
    original = original_word.lower()
    abbreviation = abbreviated_word.lower()
    words = original.split()
    initials = ''.join(word[0] for word in words if word)
    # Check if abbreviation matches initials
    if initials == abbreviation and random.random() < chance:
        return '.'.join(initials) + '.'
    # Check if original starts with abbreviation
    elif words[0].startswith(abbreviation) and words[0] != abbreviation and random.random() < chance:
        return abbreviation + '.'
    # Check if original contains abbreviation
    elif original.find(abbreviation) != -1 and words[0] != abbreviation and random.random() < chance:
        return abbreviation + '.'
    # Check if original has all abbreviation characters in order
    elif all(char in iter(original) for char in abbreviation) and random.random() < chance:
        if original.replace('-', '') == abbreviation: # auto-greaser -> autogreaser
            return abbreviation
        if abbreviation in words:
            return abbreviation
        return abbreviation + '.'
    return abbreviation

# Missing spaces in a sentence
def omit_space(sentence):
    """ Randomly omits a space from the given sentence. """
    space_idx = [idx for idx, char in enumerate(sentence) if char == ' ']
    if not space_idx: # No spaces to omit
        return sentence
    remove_idx = random.choice(space_idx)
    return sentence[:remove_idx] + sentence[remove_idx+1:]

# Extra space in a word
def add_space(word):
    """ Randomly adds a space within a word. """
    if len(word) < 3:
        return word  # Not enough characters to add a space
    index = random.randint(1, len(word) - 1)  # Ensure space is not at the beginning
    return word[:index] + ' ' + word[index:]

# Swap adjacent letters in a word
def swap_adjacent(word):
    """ Randomly swaps two adjacent letters in a given word. """
    if len(word) < 3: # Not enough letters to swap
        return word
    index = random.randint(1, len(word) - 2)
    return word[:index] + word[index + 1] + word[index] + word[index + 2:]

# Missing letter in a word
def omit_letter(word):
    """ Randomly omits one letter from a given word. """
    if len(word) < 3: # Do not omit from short words
        return word
    index = random.randint(1, len(word) - 1)
    return word[:index] + word[index + 1:]

# Double up a letter in a word
def double_letter(word):
    """ Randomly doubles one letter in a given word. """
    if len(word) < 1: # Not a word
        return word
    index = random.randint(0, len(word) - 1)
    return word[:index + 1] + word[index] + word[index + 1:]

# Replace a letter in a word with an adjacent letter (keyboard)
def adjacent_key(word):
    """ Randomly replaces a letter in a given word with an adjacent letter. """
    if len(word) < 2: # Cannot replace only letter
        return word
    index = random.randint(1, len(word) - 1)
    letter = word[index]
    if letter in KEYBOARD_DICT:
        replacement = random.choice(KEYBOARD_DICT[letter])
        return word[:index] + replacement  + word[index + 1:]
    return word

# Add adjacent letter before or after a letter in a word
def adjacent_add(word):
    """ Randomly adds an adjacent letter before or after a letter in a given word. """
    if len(word) < 2: # Not a word
        return word
    index = random.randint(1, len(word) - 1)
    letter = word[index]
    if letter in KEYBOARD_DICT:
        addition = random.choice(KEYBOARD_DICT[letter])
        if random.random() < 0.5:
            return word[:index] + addition + word[index:]           # Add before
        else:
            return word[:index + 1] + addition + word[index + 1:]   # Add after
    return word

# Replace word with its homophone
def replace_homophone(word):
    """ Replace a word with one of its homophones, if available. """
    word = word.lower()
    if word not in CMU_DICT:
        return word  # No pronunciation found
    word_pron = CMU_DICT[word][0]
    # Find homophones with the same pronunciation
    homophones = [w for w, pron in CMU_DICT.items() if pron[0] == word_pron and w != word]
    # Filter homophones with Levenshtein distance <= 1
    homophones = [w for w in homophones if levenshtein_distance(word, w) <= 1]
    if homophones: # Homophones found
        return random.choice(homophones)
    return word # No homophones found

# Introduce different typos in a sentence (default probability=0.08)
def rule_introduce_typos(sentence, chance=0.05, max_typos=3):
    """ Introduce typos in a sentence with a given probability. """
    typo_funcs = [add_space, swap_adjacent, omit_letter, double_letter, adjacent_key, adjacent_add, replace_homophone]
    typo_probs = [8, 16, 16, 17, 13, 16, 14]  # Probabilities for each typo function
    
    if random.random() < chance:
        sentence = omit_space(sentence)
    
    words = sentence.split()
    # Randomly select up to 3 words for chance to introduce typos
    typos = random.sample(range(len(words)), min(len(words), max_typos))
    for i in typos:
        if random.random() < chance:
            word = words[i]
            typo_func = random.choices(typo_funcs, weights=typo_probs, k=1)[0]
            words[i] = typo_func(word)

    return ' '.join(words)

# Introduce typos in a sentence using OpenAI GPT-4
def llm_introduce_typos(openai, sentence):
    # Chance for no typos
    if random.random() < 0.15:
        return sentence
    prompt = (
        f"Introduce a few typos into the following sentence to make it look like it was written by a human. "
        f"Use a mix of the following typo types, but avoid overdoing it. The typo types are:\n"
        f"1. Missing space between words (e.g., air conditioner -> airconditioner)\n"
        f"2. Additional space within words (e.g., permalube -> perma lube)\n"
        f"3. Swapped adjacent characters (e.g., crack -> carck)\n"
        f"4. Missing characters in a word (e.g., crack -> crak)\n"
        f"5. Double-up characters in a word (e.g., crack -> craack)\n"
        f"6. Incorrect character in a word (due to keys proximity) (e.g., crack -> xrack)\n"
        f"7. Extra characters in a word (due to keys proximity) (e.g., crack -> cracvk)\n"
        f"8. Incorrect spelling (homophones) (e.g., motor -> moter)\n\n"
        f"Here is the sentence to modify: '{sentence}'"
        f"Return the modified sentence and nothing else."
    )
    response = openai.chat.completions.create(
        model="gpt-4o-mini",
        messages=[{"role": "system", "content": "You are an expert in adding realistic typos to sentences."},
                  {"role": "user", "content": prompt}],
        temperature=0.9,
        top_p=0.9,
        n=1
    )
    if response.choices[0].message.content.startswith("'") and response.choices[0].message.content.endswith("'"):
        return response.choices[0].message.content[1:-1]
    return response.choices[0].message.content

# Humanise a MWO sentence
def humanise_sentence(sentence, llm=False):
    """ Humanise a sentence by introducing contractions, abbreviations, and typos. """
    sentence = introduce_contractions(sentence)
    sentence = introduce_abbreviations(sentence)
    if llm:
        sentence = llm_introduce_typos(llm, sentence)
    else:
        sentence = rule_introduce_typos(sentence)
    return sentence

if __name__ == '__main__':
    # Load environment variables
    load_dotenv()
    api_key = os.getenv("API_KEY")
    client = OpenAI(api_key=api_key)
    
    # Initialise global variables (dictionaries)
    current_dir = os.path.dirname(os.path.abspath("__file__"))
    main_dir = os.path.join(current_dir, '..')
    initialise_globals(main_dir)
    
    # Use humanise_sentence function
    sentence = "The air conditioner was broken."
    humanised_sentence = humanise_sentence(sentence)
    print(f"Original: {sentence}")
    print(f"Humanised: {humanised_sentence}")
    