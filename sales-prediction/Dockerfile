FROM node:14-alpine as builder

WORKDIR /

COPY ./sales-frontend/ .
RUN npm install && npm run build

FROM python:3.8-slim

WORKDIR /

RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && apt-get purge -y --auto-remove \
    && rm -rf /var/lib/apt/lists/*

COPY sales-backend/requirements.txt /requirements.txt
COPY sales-backend/requirements-base.txt /requirements-base.txt
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --compile --no-cache-dir -r /requirements.txt \
    && rm -rf /root/.cache/pip

COPY sales-backend/sentiment /sentiment

COPY sales-backend/trigger_download.py /trigger_download.py
RUN python /trigger_download.py

RUN mkdir /sales/static
COPY --from=builder /dist/ /sales/static/

CMD ["uvicorn", "sales.main:app", "--host", "0.0.0.0", "--port", "8000"]
