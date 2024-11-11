// src/routes/movies.js
const Router = require('koa-router');
const { Ingredient } = require('../models'); // Assurez-vous d'importer correctement votre modèle Ingredient

const router = new Router();

router.post('/', async (ctx) => {
  try {
    console.log('Solicitud recibida para la creación de un ingrediente');
    console.log('Datos de la solicitud:', ctx.request.body);

    const ingredient = await Ingredient.create(ctx.request.body);
    console.log('Ingrediente creado con éxito:', ingredient);

    ctx.status = 201;
    ctx.body = ingredient;
  } catch (error) {
    console.error('Error al crear el ingrediente:', error);
    ctx.status = 500;
    ctx.body = { error: 'Ocurrió un error al crear el ingrediente' };
  }
});

router.delete('/:id', async (ctx) => {
  try {
    const ingredient = await Ingredient.findByPk(ctx.params.id);
    console.log(ingredient)
    if (ingredient) {
      await ingredient.destroy();
      ctx.status = 204; // Sin contenido
    } else {
      ctx.body = { error: 'Ingrediente no encontrado' };
      ctx.status = 404; // No encontrado
    }
  } catch (error) {
    console.error('Error al eliminar el ingrediente:', error);
    ctx.body = { error: 'Ocurrió un error al eliminar el ingrediente' };
    ctx.status = 500; // Error del servidor
  }
});

// GET para obtener todos los ingredientes
router.get('/', async (ctx) => {
  try {
    const ingredients = await Ingredient.findAll(); // Consulta todos los ingredientes

    ctx.status = 200; // OK
    ctx.body = ingredients;
  } catch (error) {
    console.error('Error al obtener los ingredientes:', error);
    ctx.status = 500; // Error del servidor
    ctx.body = { error: 'Ocurrió un error al obtener los ingredientes' };
  }
});

router.patch('/:id', async (ctx) => {
  try {
    const { id } = ctx.params;
    const { body } = ctx.request;

    const ingredient = await Ingredient.findByPk(id);

    if (!ingredient) {
      ctx.status = 404; // No encontrado
      ctx.body = { error: 'Ingrediente no encontrado' };
      return;
    }

    await ingredient.update(body);

    ctx.status = 200; // OK
    ctx.body = ingredient;
  } catch (error) {
    console.error('Error al actualizar el ingrediente:', error);
    ctx.status = 500; // Error del servidor
    ctx.body = { error: 'Ocurrió un error al actualizar el ingrediente' };
  }
});

// GET para obtener ingredientes por email
router.get('/user/:email', async (ctx) => {
  const email = ctx.params.email;

  try {
    // Consultar ingredientes por email
    const ingredients = await Ingredient.findAll({
      where: {
        owner: email
      }
    });

    ctx.status = 200; // OK
    ctx.body = ingredients;
  } catch (error) {
    console.error('Error al obtener ingredientes por email:', error);
    ctx.status = 500; // Error del servidor
    ctx.body = { error: 'Ocurrió un error al obtener los ingredientes por email' };
  }
});

// GET para obtener ingredientes por ID
router.get('/:id', async (ctx) => {
  const id = ctx.params.id;

  try {
    // Consultar ingrediente por ID
    const ingredient = await Ingredient.findByPk(id);

    if (!ingredient) {
      ctx.status = 404; // No encontrado
      ctx.body = { error: 'Ingrediente no encontrado' };
      return;
    }

    ctx.status = 200; // OK
    ctx.body = ingredient;
  } catch (error) {
    console.error('Error al obtener ingrediente por ID:', error);
    ctx.status = 500; // Error del servidor
    ctx.body = { error: 'Ocurrió un error al obtener el ingrediente por ID' };
  }
});

module.exports = router;
