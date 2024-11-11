const Koa = require('koa');
const Logger = require('koa-logger');
const cors = require('@koa/cors');
const bodyParser = require('koa-bodyparser'); // Importer koa-bodyparser
const router = require('./routes'); // Importer le routeur principal

const app = new Koa();

app.use(Logger());
app.use(cors()); // Utiliser le middleware CORS
app.use(bodyParser()); // Utiliser koa-bodyparser pour parser le corps des requêtes

// Middleware pour logger toutes les requêtes
app.use(async (ctx, next) => {
  console.log(`Received ${ctx.method} request for ${ctx.url}`);
  await next();
});

app.use(router.routes()); // Utiliser le routeur principal

const PORT = process.env.PORT || 3000;
app.listen(PORT, () => {
  console.log(`Server is running on http://localhost:${PORT}`);
});
