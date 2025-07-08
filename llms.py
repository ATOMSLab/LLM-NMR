import os
from dotenv import load_dotenv
from anthropic import Anthropic
from openai import OpenAI
import google.generativeai as genai
from prompts import system_prompt


# loading variables from .env file
load_dotenv()
OAI_KEY=os.environ["OPEN_AI_KEY"]
GEMINI_KEY=os.environ["GEMINI_AI_KEY"]
ANTHROPIC_KEY=os.environ["ANTHROPIC_API_KEY"]

def call_openAI(prompt, model="gpt-4o", temperature=1):
    try:
        client = OpenAI(api_key=OAI_KEY)
        completion = client.chat.completions.create(
            model=model,
            temperature=temperature,
            messages=[ {"role": "system", "content": system_prompt},{"role": "user", "content": prompt}]
        )
        txt = completion.choices[0].message.content
        return txt,txt.split("### Start answer ###")[1].split("### End answer ###")[0].replace('\n', '').strip().lower()
    except Exception as e:
        print("[ERROR] - openAI function:\n", str(e))
        raise


def call_claude_sonnet(prompt,  model="claude-3-5-sonnet-20241022",temperature=1):
    try:
        client = Anthropic(api_key=ANTHROPIC_KEY)
        response = client.messages.create(
            model=model,
            max_tokens=1000,
            temperature=temperature,
            system=system_prompt,
            messages=[{"role": "user", "content": prompt}]
        )
        txt = response.content[0].text
        return txt,txt.split("### Start answer ###")[1].split("### End answer ###")[0].replace('\n', '').strip().lower()
    except Exception as e:
        print("[ERROR] - claude function:\n", str(e))
        raise


def call_gemini(prompt, model="gemini-2.0-flash-exp", temperature=1):
    try:
        genai.configure(api_key=GEMINI_KEY)
        model = genai.GenerativeModel(model, system_instruction=system_prompt)
        response = model.generate_content(
            contents=prompt,
            request_options={"timeout": 65000}
        )
        txt = response.text
        return txt.split("### Start answer ###")[1].split("### End answer ###")[0].replace('\n', '').strip().lower()
    except Exception as e:
        print("[ERROR] - gemini function:\n", str(e))
        raise

  