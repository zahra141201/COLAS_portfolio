const Router = require('koa-router');
const { User } = require('../models');

const router = new Router();

// Route pour créer un utilisateur
router.post('/', async (ctx) => {
  try {
    console.log('Solicitud POST recibida en /users');
    console.log('Datos recibidos:', ctx.request.body);

    const user = await User.create(ctx.request.body);
    console.log('Usuario creado:', user);

    ctx.body = user;
    ctx.status = 201;
  } catch (error) {
    console.error('Error al crear usuario:', error);
    ctx.body = error;
    ctx.status = 400;
  }
});

// Route pour obtenir tous les utilisateurs
router.get('/', async (ctx) => {
  try {
    const users = await User.findAll();
    ctx.body = users;
    ctx.status = 200;
  } catch (error) {
    console.error('Error en la consulta Sequelize:', error);
    ctx.body = error;
    ctx.status = 500;
  }
});

router.get('/:email', async (ctx) => {
  try {
    const { email } = ctx.params;
    console.log('Recherche d\'un utilisateur avec l\'adresse e-mail :', email);

    // Ajouter des logs pour vérifier ce qui est renvoyé par la base de données
    const user = await User.findOne({ where: { email } });

    console.log('Résultat de la requête de la base de données :', user);

    if (!user) {
      console.log('Aucun utilisateur trouvé avec l\'adresse e-mail :', email);
      ctx.status = 404;
      ctx.body = { error: 'No user found with the provided email address' };
    } else {
      console.log('Utilisateur trouvé avec succès :', user);
      ctx.status = 200;
      ctx.body = user;
    }
  } catch (error) {
    console.error('Une erreur est survenue lors de la recherche de l\'utilisateur par e-mail :', error);
    ctx.status = 500;
    ctx.body = { error: 'Ocurrió un error al buscar el usuario por correo electrónico' };
  }
});

// Route pour supprimer un utilisateur par email
router.delete('/deleteByEmail/:email', async (ctx) => {
  try {
    const { email } = ctx.params;
    const user = await User.findOne({ where: { email } });

    if (user) {
      await user.destroy();
      ctx.status = 204;
    } else {
      ctx.status = 404;
      ctx.body = { error: 'User not found' };
    }
  } catch (error) {
    console.log(error);
    ctx.status = 500;
    ctx.body = { error: 'An error occurred while deleting the user' };
  }
});

// Route pour mettre à jour un utilisateur par ID
router.patch('/:id', async (ctx) => {
  try {
    const { id } = ctx.params;
    const { body } = ctx.request;
    const user = await User.findByPk(id);
    if (!user) {
      ctx.status = 404;
      ctx.body = { error: 'User not found' };
      return;
    }
    await user.update(body);
    ctx.status = 200;
    ctx.body = user;
  } catch (error) {
    ctx.status = 500;
    ctx.body = { error: 'An error occurred while updating the user' };
  }
});

module.exports = router;
