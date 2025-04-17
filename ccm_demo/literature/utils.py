import os
import copy

from PIL import Image
import pytesseract
import layoutparser as lp
from layoutparser.models import Detectron2LayoutModel, detectron2

import requests

import pymupdf
from pyarrow.lib import table_to_blocks

from transformers import Qwen2_5_VLForConditionalGeneration, AutoTokenizer, AutoProcessor
from qwen_vl_utils import process_vision_info

config=os.path.abspath(os.path.join(os.path.dirname(__file__),"models/config.yaml"))
model=os.path.abspath(os.path.join(os.path.dirname(__file__),"models/model_final.pth"))

#TODO this needs to into a config for initiation
lp_model = lp.Detectron2LayoutModel(config, model,
                                    label_map={0: "Text", 1: "Title", 2: "List", 3: "Table", 4: "Figure"},
                                 extra_config=["MODEL.ROI_HEADS.SCORE_THRESH_TEST", 0.8],
                                 )

vl_model = Qwen2_5_VLForConditionalGeneration.from_pretrained(
    "Qwen/Qwen2.5-VL-7B-Instruct", torch_dtype="auto", device_map="auto", cache_dir="models",
)

processor = AutoProcessor.from_pretrained("Qwen/Qwen2.5-VL-7B-Instruct")


#TODO add tables test a paper with wester/northern blots
#TODO alternate method to include paper abstract
def process_pdf(pdf, lp_model=lp_model, vl_model=vl_model, processor=processor, zoomx=2, device="cuda", max_tokens=400):
    figure_messages = [{"role": "system", "content": [{"type": "text",
                                                       "text": """You are an expert biologist who is responsible for reading and interpreting scientific figures. For a given figure from a scientific paper
         interpret the figure. Do not provide comments on whether the figure is well done or not. Do not provide extra text on describing that you are looking at figure from a scientific figure. 
         Whenever possible very briefly describe each sections of the figure and then give an overall conclusion about what the figure tells us. 
         """}]},
                       {"role": "user", "content": [{"type": "image", "image": None, }], }]

    table_message = [{"role": "system", "content": [{"type": "text",
                                                       "text": """You are an expert biologist who is responsible for reading and interpreting scientific tables. For a given table from a scientific paper
             interpret the table. Do not provide comments on whether the table is well done or not. Do not provide extra text on describing that you are looking at table from a scientific publication. 
             Give an overall conclusion about what the tables tells us. 
             """}]},
                       {"role": "user", "content": [{"type": "image", "image": None, }], }]

    doc = pymupdf.open(pdf)
    zoom_x = zoomx  # horizontal zoom
    zoom_y = zoomx  # vertical zoom
    mat = pymupdf.Matrix(zoom_x, zoom_y)
    texts = []
    figures = []
    tables=[]
    for page in doc:
        pix = page.get_pixmap(matrix=mat)
        pix = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
        layout = lp_model.detect(pix)
        figure_blocks = lp.Layout([b for b in layout if b.type == 'Figure'])
        table_blocks = lp.Layout([b for b in layout if b.type == 'Figure'])
        if len(figure_blocks) > 0:
            for block in figure_blocks:
                coords = block.block
                coords = (coords.x_1, coords.y_1, coords.x_2, coords.y_2,)
                figure_img = pix.crop(coords)
                figures.append(figure_img)

        if len(table_blocks) > 0:
            for block in table_blocks:
                coords = block.block
                coords = (coords.x_1, coords.y_1, coords.x_2, coords.y_2,)
                table_img = pix.crop(coords)
                tables.append(figure_img)

        page_text = pytesseract.image_to_string(pix)
        texts.append(page_text)
    article_text = "\n".join(texts)

    figure_interpretation = []
    for figure in figures:
        figure_messages[1]["content"][0]["image"] = figure
        text = processor.apply_chat_template(figure_messages, tokenize=False, add_generation_prompt=True)
        # this is here for compatibility I will not be processing videos
        image_inputs, video_inputs = process_vision_info(figure_messages)
        inputs = processor(
            text=[text],
            images=image_inputs,
            videos=video_inputs,
            padding=True,
            return_tensors="pt",
        )
        inputs = inputs.to(device)
        generated_ids = vl_model.generate(**inputs, max_new_tokens=max_tokens)
        generated_ids_trimmed = [out_ids[len(in_ids):] for in_ids,
        out_ids in zip(inputs.input_ids, generated_ids)]
        output_text = processor.batch_decode(
            generated_ids_trimmed, skip_special_tokens=True, clean_up_tokenization_spaces=False)
        figure_interpretation.append(output_text)

    table_interpretation = []
    for table in tables:
        table_message[1]["content"][0]["image"] = table
        text = processor.apply_chat_template(table_message, tokenize=False, add_generation_prompt=True)
        # this is here for compatibility I will not be processing videos
        image_inputs, video_inputs = process_vision_info(table_message)
        inputs = processor(
            text=[text],
            images=image_inputs,
            videos=video_inputs,
            padding=True,
            return_tensors="pt",
        )
        inputs = inputs.to(device)
        generated_ids = vl_model.generate(**inputs, max_new_tokens=max_tokens)
        generated_ids_trimmed = [out_ids[len(in_ids):] for in_ids,
        out_ids in zip(inputs.input_ids, generated_ids)]
        output_text = processor.batch_decode(
            generated_ids_trimmed, skip_special_tokens=True, clean_up_tokenization_spaces=False)
        table_interpretation.append(output_text)

    return article_text, figures, figure_interpretation, table_interpretation

def cleanup_text(page_text):
    page_split=page_text.split("\n\n")
    sections=[]
    for section in page_split:
        cleaned=[]
        lines=section.split("\n")
        for line in lines:
            if len(line) < 50 :
                continue
            else:
                cleaned.append(line)
        sections.append("\n".join(cleaned))
    sections=[item for item in sections if len(item)>0]
    sections="\n\n".join(sections)
    return sections

#TODO image embedding model
def embed_images(figures, model):
    pass

#TODO text embeddign model
def embed_text(text, splitting_strategy, model):
    pass

