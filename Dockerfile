FROM python:3.12

WORKDIR /usr/src/app

COPY . .

RUN pip install --no-cache-dir -r requirements.txt

CMD streamlit run ./app.py --server.maxUploadSize 5000