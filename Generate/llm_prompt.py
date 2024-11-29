import os
import re
import sys
from openai import OpenAI
from dotenv import load_dotenv
from sentence_transformers import SentenceTransformer

# Initialise list of prompt variants
def initialise_prompts(openai, num_variants, num_examples):
    """ Initialise list of prompt variants. """

    # Base prompt for generating MWO sentences
    base_prompts = []
    while len(base_prompts) < num_variants:
        if num_examples == 1:
            base_prompt = "Generate a Maintenance Work Order (MWO) sentence describing the following equipment and undesirable event."
            keywords = ["Maintenance Work Order", "MWO", "equipment", "undesirable event", "sentence"]
        elif num_examples == 5:
            base_prompt = f"Generate {num_examples} different Maintenance Work Order (MWO) sentences describing the following equipment and undesirable event."
            keywords = [f"{num_examples}", "Maintenance Work Order", "MWO", "equipment", "undesirable event", "sentence"]
        base_prompts = paraphrase_prompt(openai, base_prompt, keywords, num_variants)
        similarity = check_similarity(base_prompt, base_prompts)
        for prompt, sim in zip(base_prompts, similarity):
            if sim > 0.9:
                base_prompts.append(prompt)
        base_prompts = list(set(base_prompts)) # Remove duplicates

    # Word verbose-limit instruction for generating MWO sentences
    limit_words = []
    while len(limit_words) < num_variants:
        instruction = "Avoid verbosity and use minimal stop words."
        keywords = ["verbosity", "stop words"]
        limit_words = paraphrase_prompt(openai, instruction, keywords, num_variants)
        similarity = check_similarity(instruction, limit_words)
        for prompt, sim in zip(limit_words, similarity):
            if sim > 0.9:
                limit_words.append(prompt)
        limit_words = list(set(limit_words)) # Remove duplicates

    # Word count-limit instruction for generating MWO sentences
    limit_count = []
    while len(limit_count) < num_variants:
        if num_examples == 1:
            instruction = "The sentence can have a maximum of 8 words."
        elif num_examples == 5:
            instruction = "Each sentence can have a maximum of 8 words."
        keywords = ["sentence", "8"]
        limit_count = paraphrase_prompt(openai, instruction, keywords, num_variants)
        similarity = check_similarity(instruction, limit_count)
        for prompt, sim in zip(limit_count, similarity):
            if sim > 0.9:
                limit_count.append(prompt)
        limit_count = list(set(limit_count)) # Remove duplicates

    # Make sure list has correct number of variants
    base_prompts = base_prompts[:num_variants]
    limit_words = limit_words[:num_variants]
    limit_count = limit_count[:num_variants]
    return (base_prompts, limit_words, limit_count)

# Check semantic similarity for the paraphrased sentences
def check_similarity(original, paraphrases):
    """ Check semantic similarity for the paraphrased sentences. """
    model = SentenceTransformer('sentence-transformers/paraphrase-MiniLM-L6-v2')
    original_embedding = model.encode(original)
    paraphrases_embeddings = model.encode(paraphrases)
    similarities = model.similarity(original_embedding, paraphrases_embeddings)
    
    # Uncomment this to print sentence and their similarity scores
    # for sentence, similarity in zip(paraphrases, similarities.tolist()[0]):
    #     print(f"{similarity:.4f} - {sentence}")
    return similarities.tolist()[0]

# Post-process LLM response of prompt paraphrases into a list of sentences
def process_prompt_response(response):
    """ Process the response from the LLM. """
    output = []
    sentences = response.split('\n')
    for sentence in sentences:
        processed = re.sub(r'^\d+\.\s*', '', sentence).strip() 
        output.append(processed)
    return output

# Get LLM to paraphrase prompts for generating more diverse responses
def paraphrase_prompt(openai, prompt, keywords=None, num_paraphrases=5):
    """ Paraphrase the prompt to generate more diverse responses. """
    # Paraphrase the prompt num_paraphrases times
    paraphrase_prompt = f"Paraphrase the following sentence {num_paraphrases} times.\n{prompt}\n"
    paraphrase_prompt += "\n" + "Do not add any new information or alter the meaning."
    if keywords:
        string_keywords = ", ".join(keywords)
        paraphrase_prompt += "Must include the following keywords: " + string_keywords
    response = openai.chat.completions.create(
                    model="gpt-4o-mini",
                    messages=[
                            {"role": "system", "content": "You are a sentence paraphraser."},
                            {"role": "user", "content": paraphrase_prompt},
                        ],
                    top_p=0.9,
                    temperature=0.9,
                    n=1
                )
    output = process_prompt_response(response.choices[0].message.content)
    return output

if __name__ == "__main__":
    # Set OpenAI API key
    load_dotenv()
    api_key = os.getenv("API_KEY")
    client = OpenAI(api_key=api_key)
    
    # Test Paraphrase Instruction Prompt
    instruction = "The sentence can have a maximum of 8 words."
    instructions_dummy = [
        "The sentence can have a maximum of 8 words.",
        "Each sentence can have a maximum of 8 words.",
        "The sentence should have a maximum of 8 words.",
        "Each sentence should have a maximum of 8 words.",
        "The sentence must have a maximum of 8 words."
    ]
    similarity = check_similarity(instruction, instructions_dummy)
    for prompt, sim in zip(instructions_dummy, similarity):
        print(f"{sim:.4f} - {prompt}")