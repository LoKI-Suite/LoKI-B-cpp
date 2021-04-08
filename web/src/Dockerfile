FROM node:alpine
RUN apk add python3
WORKDIR /app
COPY . .
CMD ["python3", "-m", "http.server", "3000"]
