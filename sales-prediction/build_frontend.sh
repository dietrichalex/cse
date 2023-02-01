#!/usr/bin/env bash

cd ./sales-frontend/ || exit 1
npm run build
cd ..

rm -r ./sales-backend/sales/static/* || mkdir ./sales-backend/sales/static
cp -r ./sales-frontend/dist/* ./sales-backend/sales/static/
