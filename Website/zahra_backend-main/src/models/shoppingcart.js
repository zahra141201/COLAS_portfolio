'use strict';
const { Model } = require('sequelize');

module.exports = (sequelize, DataTypes) => {
  class ShoppingCart extends Model {
    static associate(models) {
      this.belongsTo(models.User, {
        foreignKey: 'email_cliente',
        onDelete: 'CASCADE', 
        onUpdate: 'CASCADE',// Supprime les paniers associés lorsque l'utilisateur est supprimé
      });
      this.belongsTo(models.Ingredient, {
        foreignKey: 'id_ingredient',
        onDelete: 'CASCADE', 
        onUpdate: 'CASCADE',// Supprime les paniers associés lorsque l'ingrédient est supprimé
      });
    }
  }

  ShoppingCart.init({
    email_cliente: DataTypes.STRING,
    id_ingredient: DataTypes.INTEGER,
    estado: DataTypes.BOOLEAN
  }, {
    sequelize,
    modelName: 'ShoppingCart',
  });

  return ShoppingCart;
};
